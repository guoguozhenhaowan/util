#include "onlinebwa.h"

OnlineBWA::OnlineBWA(){
    mIndex = 0;
    mCopyComment = false;
    mMemOpt = mem_opt_init();
    mMemOpt->flag |= MEM_F_SOFTCLIP;
}

OnlineBWA::~OnlineBWA(){
    if(mIndex){
        bwa_idx_destroy(mIndex);
    }
    if(mMemOpt){
        free(mMemOpt);
    }
}

bam_hdr_t* OnlineBWA::getBamHeader(){
    if(!mIndex){
        return NULL;
    }
    bntseq_t* bns = mIndex->bns;
    std::stringstream ss;
    for(int i = 0; i < bns->n_seqs; ++i){
        ss << "@SQ\tSN:" << bns->anns[i].name << "\tLN:" << bns->anns[i].len << "\n";
    }
    kstring_t str;
    bam_hdr_t *hhh;
    str.l = str.m = 0; str.s = 0;
    std::istringstream iss(ss.str());
    std::string line;
    while(std::getline(iss, line, '\n')){
        if(line.length() == 0 || line.at(0) != '@'){
            break;
        }
        kputsn(line.c_str(), line.length(), &str);
        kputc('\n', &str);
    }
    if(str.l == 0){
        kputsn("", 0, &str);
    }
    hhh = sam_hdr_parse(str.l, str.s);
    hhh->l_text = str.l; hhh->text = str.s; // hhh->text needs to be freed
    return hhh;
}

bwt_t* OnlineBWA::pac2bwt(const uint8_t* pac, int bwt_seq_lenr){
    uint64_t i = 0;
    bwt_t* bwt = (bwt_t*)std::calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwt_seq_lenr;
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;
    std::memset(bwt->L2, 0, 5 * 4);
    ubyte_t* buf = (ubyte_t*)std::calloc(bwt->seq_len + 1, 1);
    for(i = 0; i < bwt->seq_len; ++i){
        buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
        ++bwt->L2[1+buf[i]];
    }
    for(i = 2; i <= 4; ++i){
        bwt->L2[i] += bwt->L2[i-1];
    }
    bwt->primary = is_bwt(buf, bwt->seq_len);
    bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
    for(i = 0; i < bwt->seq_len; ++i){
        bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
    }
    free(buf);
    return bwt;
}

void OnlineBWA::addAnns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset){
    ann->offset = offset;
    ann->name = (char*)std::malloc(name.length() + 1);
    std::strncpy(ann->name, name.c_str(), name.length() + 1);
    ann->anno = (char*)std::malloc(7);
    std::strncpy(ann->anno, "(null)\0", 7);
    ann->len = seq.length();
    ann->n_ambs = 0;
    ann->gi = 0;
    ann->is_alt = 0;
}

uint8_t* OnlineBWA::addSeq2pac(const kseq_t* seq, bntseq_t* bns, uint8_t* pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q){
    bntann1_t* p;
    int lasts;
    if(bns->n_seqs == *m_seqs){
        *m_seqs <<= 1;
        bns->anns = (bntann1_t*)std::realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
    }
    p = bns->anns + bns->n_seqs;
    p->name = strdup((char*)seq->name.s);
    p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
    p->gi = 0; p->len = seq->seq.l;
    p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
    p->n_ambs = 0;
    for(size_t i = lasts = 0; i < seq->seq.l; ++i) {
        int c = nst_nt4_table[(int)seq->seq.s[i]];
        if(c >= 4) { // N
            if (lasts == seq->seq.s[i]) { // contiguous N
                ++(*q)->len;
            }else{
                if(bns->n_holes == *m_holes){
                    (*m_holes) <<= 1;
                    bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                }
                *q = bns->ambs + bns->n_holes;
                (*q)->len = 1;
                (*q)->offset = p->offset + i;
                (*q)->amb = seq->seq.s[i];
                ++p->n_ambs;
                ++bns->n_holes;
            }
        }
        lasts = seq->seq.s[i];
        if (c >= 4) c = lrand48()&3;
        if (bns->l_pac == *m_pac) { // double the pac size
        *m_pac <<= 1;
        pac = (uint8_t*)realloc(pac, *m_pac/4);
        memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
        _set_pac(pac, bns->l_pac, c);
        ++bns->l_pac;
        }
    }
    ++bns->n_seqs;
    return pac;
}

uint8_t* OnlineBWA::makePac(std::vector<UnalignedSeq> v, bool addReverse){
    bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    uint8_t *pac = 0;
    int32_t m_seqs, m_holes;
    int64_t m_pac, l;
    bntamb1_t *q;
    bns->seed = 11; // fixed seed for random generator
    m_seqs = m_holes = 8; m_pac = 0x10000;
    bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
    bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
    pac = (uint8_t*) calloc(m_pac/4, 1);
    q = bns->ambs;
    // move through the unaligned sequences
    for(size_t k = 0; k < v.size(); ++k){
        // make the ref name kstring
        kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
        name->l = v[k].mName.length() + 1;
        name->m = v[k].mName.length() + 3;
        name->s = (char*)std::calloc(name->m, sizeof(char));
        memcpy(name->s, v[k].mName.c_str(), v[k].mName.length()+1);
        // make the sequence kstring
        kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
        t->l = v[k].mSeq.length();
        t->m = v[k].mSeq.length() + 2;
        //t->s = (char*)calloc(v[k].Seq.length(), sizeof(char));
        t->s = (char*)malloc(t->m);
        std::memcpy(t->s, v[k].mSeq.c_str(), v[k].mSeq.length());
        // put into a kstring
        kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));
        ks->seq = *t;
        ks->name = *name;
        // make the forward only pac
        pac = addSeq2pac(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
        // clear it out
        free(name->s);
        free(name);
        free(t->s);
        free(t);
        free(ks);
    }
    // if addReverse needed
    if(addReverse){
        //add the reverse complemented sequence
        m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
        pac = (uint8_t*)realloc(pac, m_pac/4);
        memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
        for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
       _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
    }
    bns_destroy(bns);
    return pac;
}

void OnlineBWA::writePacToFile(const std::string& file) const{
    FILE* fp;
    std::string fname = file + ".pac";
    fp = xopen(fname.c_str(), "wb");
    ubyte_t ct;
    err_fwrite(mIndex->pac, 1, (mIndex->bns->l_pac>>2) + ((mIndex->bns->l_pac&3) == 0? 0 : 1), fp);
    if(mIndex->bns->l_pac % 4 == 0){
        ct = 0;
        err_fwrite(&ct, 1, 1, fp);
    }
    ct = mIndex->bns->l_pac % 4;
    err_fwrite(&ct, 1, 1, fp);
    err_fflush(fp);
    err_fclose(fp);
}

void OnlineBWA::alignSeq(const std::string& seq, const std::string& name, std::vector<bam1_t*>& result){
    if(!mIndex){
        return;
    }
    result.clear();
    mem_alnreg_v ar = mem_align1(mMemOpt, mIndex->bwt, mIndex->bns, mIndex->pac, seq.length(), seq.data());
    for(size_t i = 0; i < ar.n; ++i){
        mem_aln_t a = mem_reg2aln(mMemOpt, mIndex->bns, mIndex->pac, seq.length(), seq.c_str(), &ar.a[i]);
        bam1_t* b = bam_init1();
        b->core.tid = a.rid;
        b->core.pos = a.pos;
        b->core.qual = a.mapq;
        b->core.flag = a.flag;
        b->core.n_cigar = a.n_cigar;
        b->core.mtid = -1;
        b->core.mpos = -1;
        b->core.isize = 0;
        if(a.is_rev){
            b->core.flag |= BAM_FREVERSE;
        }
        b->core.l_qname = name.length() + 1;
        b->core.l_qseq = seq.length();
        b->l_data = b->core.l_qname + (a.n_cigar << 2) + ((b->core.l_qseq + 1) >> 1) + (b->core.l_qseq);
        b->data = (uint8_t*)std::malloc(b->l_data);
        std::memcpy(b->data, name.c_str(), name.length() + 1);
        std::memcpy(b->data + b->core.l_qname, (uint8_t*)a.cigar, a.n_cigar << 2);
        uint8_t* mbases = b->data + b->core.l_qname + (b->core.n_cigar << 2);
        std::string fseq = seq;
        if(a.is_rev){
            fseq = util::reverseComplete(seq);
        }
        for(uint32_t j = 0; j < fseq.length(); ++j){
            uint8_t base = 15;
            switch(fseq[j]){
                case 'A':
                    base = 1;
                    break;
                case 'C':
                    base = 2;
                    break;
                case 'G':
                    base = 4;
                    break;
                case 'T':
                    base = 8;
                    break;
                default:
                    base = 15;
            }
            mbases[j >> 1] &= ~(0xF << ((~j & 1) << 2));
            mbases[j >> 1] |= base << ((~j & 1) << 2);
        }
        uint8_t* quals = bam_get_qual(b);
        quals[0] = 0xff;
        size_t arn = ar.n;
        bam_aux_append(b, "NA", 'i', 4, (uint8_t*)(&arn));
        uint32_t irec;
        irec = a.NM;
        bam_aux_append(b, "NM", 'i', 4, (uint8_t*)(&irec));
        if(a.XA){
            bam_aux_append(b, "XA", 'Z', std::strlen(a.XA), (uint8_t*)a.XA);
        }
        irec = a.score;
        bam_aux_append(b, "AS", 'i', 4, (uint8_t*)(&irec));
        result.push_back(b);
        free(a.cigar);
    }
    free(ar.a);
}

void OnlineBWA::alignSeq(const UnalignedSeq& us, std::vector<bam1_t*>& result){
    alignSeq(us.mSeq, us.mName, result);
    if(mCopyComment){
        for(uint32_t i = 0; i < result.size(); ++i){
            bam_aux_append(result[i], "BC", 'Z', us.mComment.length(), (uint8_t*)us.mComment.c_str());
        }
    }
}

void OnlineBWA::alignSeq(const kseq_t* seq, std::vector<bam1_t*>& result){
    alignSeq(seq->seq.s, seq->name.s, result);
    if(mCopyComment){
        for(uint32_t i = 0; i < result.size(); ++i){
            bam_aux_append(result[i], "BC", 'Z', seq->comment.l, (uint8_t*)seq->comment.s);
        }
    }
}

void OnlineBWA::constructIndex(const std::vector<UnalignedSeq>& uv){
    if(!uv.size()){
        return;
    }
    for(auto& e: uv){
        if(e.mName.empty() || e.mSeq.empty()){
            util::error_exit("nama and sequence must all be unempty to construct a index");
        }
    }
    if(mIndex){
        bwa_idx_destroy(mIndex);
        mIndex = 0;
    }
    // allocate memory for index
    mIndex = (bwaidx_t*)std::calloc(1, sizeof(bwaidx_t));
    // construct the forward-only pac
    uint8_t* fwd_pac = makePac(uv, false);
    // construct the forward-reverse pac
    uint8_t* pac = makePac(uv, true);
    size_t tlen = 0;
    for(auto& e: uv){
        tlen += e.mSeq.length();
    }
    // make the bwt
    bwt_t* bwt = pac2bwt(pac, tlen * 2);
    bwt_bwtupdate_core(bwt);
    free(pac);
    // construct sa from bwt and occ. add it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);
    bntseq_t* bns = (bntseq_t*)std::calloc(1, sizeof(bntseq_t));
    bns->l_pac = tlen;
    bns->n_seqs = uv.size();
    bns->seed = 11;
    bns->n_holes = 0;
    // make the anns
    bns->anns = (bntann1_t*)std::calloc(uv.size(), sizeof(bntann1_t));
    size_t offset = 0;
    for(size_t k = 0; k < uv.size(); ++k){
        addAnns(uv[k].mName, uv[k].mSeq, &bns->anns[k], offset);
        offset += uv[k].mSeq.length();
    }
    // ambs is 'holes', like N bases
    bns->ambs = 0;
    // make the in-memory idx struct
    mIndex->bwt = bwt;
    mIndex->bns = bns;
    mIndex->pac = fwd_pac;
    return;
}

void OnlineBWA::loadIndex(const std::string& file){
    bwaidx_t* newIndex = bwa_idx_load(file.c_str(), BWA_IDX_ALL);
    if(!newIndex){
        util::error_exit("error loading index");
    }
    if(newIndex){
        bwa_idx_destroy(mIndex);
    }
    mIndex = newIndex;
}

void OnlineBWA::writeIndex(const std::string& file){
    if(!mIndex){
        return;
    }
    std::string bwt_file = file + ".bwt";
    std::string sa_file = file + ".sa";
    bwt_dump_bwt(bwt_file.c_str(), mIndex->bwt);
    bwt_dump_sa(sa_file.c_str(), mIndex->bwt);
    bns_dump(mIndex->bns, file.c_str());
    writePacToFile(file);
}

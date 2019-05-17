#include "bamutil.h"

void bamutil::dump(bam1_t* b){
    std::cerr << "R:     " << b->core.tid << ":" << b->core.pos << std::endl;
    std::cerr << "M:     " << b->core.mtid << ":" <<  b->core.mpos << std::endl;
    std::cerr << "TLEN:  " << b->core.isize << std::endl;
    std::cerr << "QName: " << getQName(b) << std::endl;
    std::cerr << "Cigar: " << getCigar(b) << std::endl;
    std::cerr << "Seq:   " << getSeq(b) << std::endl;
    std::cerr << "Qual:  " << getQual(b) << std::endl;
}

std::string bamutil::getQName(bam1_t* b){
    return bam_get_qname(b);
}

std::string bamutil::getQual(bam1_t* b){
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    std::string qs(len, '\0');
    for(int i = 0; i < len; ++i){
        qs[i] = (char)(data[i] + 33);
    }
    return qs;
}

int bamutil::getED(bam1_t* b){
    const char tagNM[2] = {'N', 'M'};
    uint8_t* dataNM = bam_aux_get(b, tagNM);
    if(!dataNM){
        return 0;
    }
    return bam_aux2i(dataNM);
}

std::string bamutil::getSeq(bam1_t* b){
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    std::string seq(len, '\0');
    for(int i = 0; i < len; ++i){
        seq[i] = fourbits2base(bam_seqi(data, i));
    }
    return seq;
}

std::string bamutil::getBarcode(bam1_t* b){
    const char tagBC[2] = {'O', 'X'};
    uint8_t* dataBC = bam_aux_get(b, tagBC);
    if(!dataBC){
        return "";
    }
    return bam_aux2Z(dataBC);
}

std::string bamutil::getCigar(bam1_t* b){
    std::stringstream ss;
    uint32_t* data = bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    for(int i = 0; i < cigarNum; ++i){
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << op << len;
    }
    return ss.str();
}

bool bamutil::isPartOf(bam1_t* part, bam1_t* whole, bool isLeft){
    if(!part || !whole){
        return false;
    }
    uint32_t* cigarPart = bam_get_cigar(part);
    int cigarNumPart = part->core.n_cigar;
    uint32_t* cigarWhole = bam_get_cigar(whole);
    int cigarNumWhole = part->core.n_cigar;
    if(cigarNumWhole < cigarNumPart){
        return false;
    }

    for(int i = 0; i < cigarNumPart; ++i){
        uint32_t valPart = cigarPart[i];
        if(!isLeft){
            valPart = cigarPart[cigarNumPart - i - 1];
        }
        char opPart = bam_cigar_opchr(valPart);
        uint32_t lenPart = bam_cigar_oplen(valPart);

        uint32_t valWhole = cigarWhole[i];
        if(!isLeft){
            valWhole = cigarWhole[cigarNumWhole - i - 1];
        }
        char opWhole = bam_cigar_opchr(valWhole);
        uint32_t lenWhole = bam_cigar_oplen(valWhole);
        
        if(opPart != opWhole || lenPart > lenWhole){
            return false;
        }
        if(lenPart < lenWhole){
            if(i == cigarNumPart - 1){
                return true;
            }else{
                if(i != cigarNumPart - 2){
                    return false;
                }else{
                    int next = i + 1;
                    uint32_t valPartNext = cigarPart[next];
                    if(!isLeft){
                        valPartNext = cigarPart[cigarNumPart - next - 1];
                    }
                    int opPartNext = bam_cigar_op(valPartNext);
                    if(opPartNext != BAM_CHARD_CLIP){
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void bamutil::dumpHeader(bam_hdr_t* hdr){
    std::cerr << hdr->n_targets << " contigs in the bam file:\n";
    int dumped = 0;
    while(dumped < hdr->n_targets){
        char* targetName = hdr->target_name[dumped];
        int targetLen = hdr->target_len[dumped];
        std::cout << targetName << ": " << targetLen << " bp\n";
        ++dumped;
    }
}

int bamutil::getRefOffset(bam1_t* b, int offset){
    const int QUERY_CONSUME[16] = {1, 1, 0, 0, 1, 0, 0, 1, 1, 0};
    const int REFERENCE_CONSUME[16] = {1, 0, 1, 1, 0, 0, 0, 1, 1, 0};
    uint32_t* data = bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    int ref = 0;
    int query = 0;
    for(int i = 0; i < cigarNum; ++i){
        uint32_t val = data[i];
        int op = bam_cigar_op(val);
        uint32_t len = bam_cigar_oplen(val);
        query += len * QUERY_CONSUME[op];
        ref += len * REFERENCE_CONSUME[op];
        if(query >= offset){
            if(op == BAM_CINS || op == BAM_CSOFT_CLIP){
                return -1;
            }else{
                return ref - REFERENCE_CONSUME[op] * (query - offset);
            }
        }
    }
    std::cerr << "wrong cigar: " << getCigar(b);
    std::cerr << " tid: " << b->core.tid << " l_qseq: " << b->core.l_qseq;
    std::cerr << " pos: "  << b->core.l_qseq << " isize: " << b->core.isize;
    std::cerr << " offset: " << offset << std::endl;
    // not found
    return -1;
}

void bamutil::getFirstMatchOffsetAndLen(bam1_t* b, int& MatchOffset, int& MatchLen){
    const int QUERY_CONSUME[16] = {1, 1, 0, 0, 1, 0, 0, 1, 1, 0};
    uint32_t* data = bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    int query = 0;
    for(int i = 0; i < cigarNum; ++i){
        uint32_t val = data[i];
        int op = bam_cigar_op(val);
        uint32_t len = bam_cigar_oplen(val);
        if(op == BAM_CMATCH){
            MatchOffset = query;
            MatchLen = len;
            return;
        }
        query += len * QUERY_CONSUME[op];
    }
    MatchOffset = -1;
    MatchLen = 0;
}

uint8_t bamutil::getMappedBases(bam1_t* b){
    uint8_t mappedBases = 0;
    uint32_t* data = bam_get_cigar(b);
    uint32_t cigarNum = b->core.n_cigar;
    for(uint32_t i = 0; i < cigarNum; ++i){
        uint32_t val = data[i];
        int op = bam_cigar_op(val);
        if(op == BAM_CMATCH){
            mappedBases += bam_cigar_oplen(val);
        }
    }
    return mappedBases;
}

void bamutil::changeQname(bam1_t* b, const std::string n){
    size_t nonq_len = b->l_data - b->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    memcpy(nonq, b->data + b->core.l_qname, nonq_len);
    free(b->data);
    b->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    memcpy(b->data, (uint8_t*)n.c_str(), n.length());
    b->data[n.length()] = '\0';
    b->l_data = b->l_data - b->core.l_qname + n.length() + 1;
    b->core.l_qname = n.length() + 1;
    memcpy(b->data + b->core.l_qname, nonq, nonq_len);
    free(nonq);
    b->m_data = b->l_data;
    b->core.l_extranul = 0;
}

void bamutil::changeSeq(bam1_t* b, int i, char ch){
    uint8_t base = bamutil::base2fourbits(ch);
    uint8_t* seq = bam_get_seq(b);
    if(i < 0 || i > b->core.l_qseq - 1){
        util::error_exit("invalid position to change!");
    }
    if(i % 2 == 1){
        seq[i/2] = (seq[i/2] & 0xF0) | base;
    }else{
        seq[i/2] = (seq[i/2] & 0x0F) | (base << 4);
    }
}

void bamutil::changeSeq(bam1_t* b, int i, uint8_t ch){
    uint8_t* seq = bam_get_seq(b);
    if(i < 0 || i > b->core.l_qseq - 1){
        util::error_exit("invalid position to change!");
    }
    if(i % 2 == 1){
        seq[i/2] = (seq[i/2] & 0xF0) | ch;
    }else{
        seq[i/2] = (seq[i/2] & 0x0F) | (ch << 4);
    }
}

void bamutil::copyQName(bam1_t* from, bam1_t* to){
    char* fromName = bam_get_qname(from);
    char* toName = bam_get_qname(to);
    int fromlen = from->core.l_qname;
    int tolen = to->core.l_qname;
    if(tolen >= fromlen){
        std::memcpy(toName, fromName, fromlen);
        to->core.l_extranul = from->core.l_extranul;
        if(fromlen != tolen){
            char* start = toName + fromlen;
            char* end = (char*)to->data + to->l_data;
            int offset = tolen - fromlen;
            for(char* p = start; p + offset < end; ++p){
                *p = *(p + offset);
            }
            to->core.l_qname = fromlen;
            to->l_data -= offset;
        }
    }else{
        uint8_t* newdata = new uint8_t[fromlen + to->l_data - tolen];
        std::memcpy(newdata, fromName, fromlen);
        std::memcpy(newdata + fromlen, to->data + tolen, to->l_data - tolen);
        std::free(to->data);
        to->data = newdata;
        to->core.l_extranul = from->core.l_extranul;
        to->core.l_qname = fromlen;
        to->l_data = fromlen + to->l_data - tolen;
    }
}

bool bamutil::isPrimary(bam1_t* b){
    if(b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY){
        return false;
    }
    return true;
}

bool bamutil::isProperPair(bam1_t* b){
    return b->core.flag & BAM_FPROPER_PAIR;
}

int bamutil::getRightRefPos(bam1_t* b){
    return bam_endpos(b);
}

uint8_t bamutil::base2fourbits(char base){
    switch(base){
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 4;
        case 'T':
            return 8;
        case 'N':
            return 15;
        default:
            return 15;
    }
}

char bamutil::fourbits2base(uint8_t val){
    switch(val){
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            return 'N';
    }
}

uint8_t bamutil::getReadQualSum(bam1_t* b){
    if(!b){
        return 0;
    }
    uint8_t sum = 0;
    uint8_t* qseq = bam_get_qual(b);
    for(int32_t i = 0; i < b->core.l_qseq; ++i){
        sum += qseq[i];
    }
    return sum;
}

std::vector<int32_t> bamutil::getMismatchPos(bam1_t* b){
    std::vector<int32_t> misp;
    std::vector<std::string> nums;
    uint8_t* mdTag = bam_aux_get(b, "MD");
    if(!mdTag){
        return misp;
    }
    std::string mdStr = bam_aux2Z(mdTag);
    std::string::size_type las, cur;
    las = mdStr.find_first_of("0123456789");
    // first pass, find all matched lengths
    while((las = mdStr.find_first_of("0123456789", las)) != std::string::npos){
        cur = mdStr.find_first_not_of("0123456789", las);
        if(cur != std::string::npos){
            nums.push_back(mdStr.substr(las, cur - las));
        }else{
            nums.push_back(mdStr.substr(las));
            break;
        }
        las = cur;
    }
    // second pass, get mismatched positions
    int32_t clen = 0;
    las = 0;
    for(uint32_t i = 0; i < nums.size(); ++i){
        cur = mdStr.find(nums[i], las);
        if(cur == 0){
            clen += std::atoi(nums[i].c_str());
            las = cur + nums[i].length();
            continue;
        }
        std::string lasStr = mdStr.substr(las, cur - las);
        if(lasStr.find("^") != std::string::npos){
            clen += lasStr.length() - 1;
        }else if(lasStr.find_first_of("ATCG") != std::string::npos){
            misp.push_back(clen);
            clen += 1;
        }
        clen += std::atoi(nums[i].c_str());
        las = cur + nums[i].length();
    }
        std::string fstr = mdStr.substr(las);
    if(!fstr.empty()){
        if(fstr.find("^") == std::string::npos && fstr.find_first_of("ATCG") != std::string::npos){
            misp.push_back(clen);
        }
    }
    for(uint32_t i = 0; i < misp.size(); ++i){
        misp[i] += b->core.pos;
    }
    return misp;
}

void bamutil::getIndelPos(const bam1_t* b, std::vector<int>& indelPosVec){
    indelPosVec.clear();
    uint32_t offset = 0;
    uint32_t* data = bam_get_cigar(b);
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        char op = bam_cigar_opchr(data[i]);
        uint32_t len = bam_cigar_oplen(data[i]);
        // skip cigar which does not consume query except BAM_CDEL
        if(op == BAM_CPAD || op == BAM_CHARD_CLIP || op == BAM_CREF_SKIP){
            continue;
        }
        // increase offset only if cigar consume query except BAM_CINS
        if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CSOFT_CLIP){
            offset += len;
            continue;
        }
        // record BAM_CINS
        if(op == BAM_CINS){
            indelPosVec.push_back(offset);
            offset += len;
            indelPosVec.push_back(offset);
            continue;
        }
        // record BAM_CDEL
        if(op == BAM_CDEL){
            indelPosVec.push_back(offset);
            continue;
        }
    }
}

bool bamutil::posIsAdjacentIndel(int pos, int range, const std::vector<int>& ivec){
    for(uint32_t i = 0; i < ivec.size(); ++i){
        if(std::abs(pos - ivec[i]) <= range){
            return true;
        }
    }
    return false;
}

int bamutil::getCompactMismatch(bam1_t* b){
    int inDelLen = 0;
    int inDelNum = 0;
    int nBaseNum = 0;
    int mismatchNum = 0;
    std::string seq = getSeq(b);
    // ambiguous N base would not be counted as mismatch
    for(uint32_t i = 0; i < seq.length(); ++i){
        if(seq[i] == 'N'){
            ++nBaseNum;
            --mismatchNum;
        }
    }
    uint32_t* cigar = bam_get_cigar(b);
    bool cigarMGot = false;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        int32_t op = bam_cigar_op(cigar[i]);
        int32_t ol = bam_cigar_oplen(cigar[i]);
        switch(op){
            case BAM_CMATCH:
                cigarMGot = true;
                break;
            case BAM_CDIFF:
                ++mismatchNum;
                break;
            case BAM_CINS: case BAM_CDEL:
                ++mismatchNum; // each indel counted as one mismatch
                inDelLen += ol;
                ++inDelNum;
                break;
            default:
                break;
        }
    }
    // adjust cigar M
    if(cigarMGot){
        uint8_t* nmTag = (uint8_t*)bam_aux_get(b, "NM");
        if(nmTag){
            mismatchNum = bam_aux2i(nmTag) - inDelLen + inDelNum - nBaseNum;
        }
    }
    return mismatchNum < 0 ? 0 : mismatchNum;
}

int bamutil::getBamReadLength(const std::string& bamFile){
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    int readLen = 0;
    int count = 0;
    while(count < 10000 && sam_read1(fp, h, b) >= 0){
        if(b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL)){
            continue;
        }
        if(b->core.n_cigar == 1){
            readLen = std::max(b->core.l_qseq, readLen);
            ++count;
        }
    }
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
    return readLen;
}

bool bamutil::readRearPartPileup(const bam_pileup1_t* p, int range){
    if(bam_is_rev(p->b)){
        if(p->qpos <= range){
            return true;
        }
    }else{
        if(p->qpos >= p->b->core.l_qseq - range){
            return true;
        }
    }
    return false;
}

float bamutil::getAverageBaseQualityAroundPileupPos(const bam_pileup1_t* p, int range){
    uint8_t* qualStr = bam_get_qual(p->b);
    float totalQual = 0.0;
    int i = 0;
    for(i = std::max(0, p->qpos - range); i <= std::min(p->b->core.l_qseq, p->qpos + range); ++i){
        totalQual += qualStr[i];
    }
    return totalQual/i;
}

std::string bamutil::getInsertionSeq(const bam_pileup1_t* p, int maxLen){
    int insLenToGet = std::min(p->indel, maxLen);
    std::string ret(insLenToGet, '\0');
    for(int i = 0; i < insLenToGet; ++i){
        ret[i] = fourbits2base(bam_seqi(bam_get_seq(p->b), p->qpos + i));
    }
    return ret;
}

std::pair<int, int> bamutil::getSoftClipLength(bam1_t* b){
    std::pair<int, int> ret = {0, 0};
    uint32_t* cigarStr = bam_get_cigar(b);
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        int op = bam_cigar_op(cigarStr[i]);
        int ol = bam_cigar_oplen(cigarStr[i]);
        if(op == BAM_CSOFT_CLIP){
            if(i == 0){
                ret.first = ol;
            }else{
                ret.second = ol;
            }
        }
    }
    return ret;
}

int bamutil::getMinDistToReadEnd(const bam_pileup1_t* p){
    std::pair<int, int> clipLen = getSoftClipLength(p->b);
    return std::min(p->qpos + 1 - clipLen.first, p->b->core.l_qseq - clipLen.second - p->qpos - 1);
}

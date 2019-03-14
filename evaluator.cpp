#include "evaluator.h"

int Evaluator::seq2int(std::string& seq, int pos, int keylen, int lastVal){
    int rlen = seq.length();
    if(lastVal >= 0){
        const int mask = (1 << (keylen * 2)) - 1;
        int key = (lastVal << 2) & mask;
        char base = seq[pos + keylen - 1];
        switch(base){
            case 'A':
                key += 0;
                break;
            case 'T':
                key += 1;
                break;
            case 'C':
                key += 2;
                break;
            case 'G':
                key += 3;
                break;
            default:
                return -1;
        }
        return key;
    }else{
        int key = 0;
        for(int i = pos; i < pos + keylen; ++i){
            key = (key << 2);
            char base = seq[i];
            switch(base){
                case 'A':
                    key += 0;
                    break;
                case 'T':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                default:
                    return -1;
            }
        }
        return key;
    }
}

std::string Evaluator::int2seq(size_t val, int seqLen){
    char bases[4] = {'A', 'T', 'C', 'G'};
    std::string ret(seqLen, 'N');
    int index = 0;
    while(index < seqLen){
        ret[seqLen - index -1] = bases[val & 0x03];
        val = (val >> 2);
        ++index;
    }
    return ret;
}

bool Evaluator::isTwoColorSystem(){
    FqReader fqr(this->r1File);
    Read* r = fqr.read();

    if(!r){
        return false;
    }
    // NEXTSEQ500, NEXTSEQ 550, NOVASEQ are two color system and with specific fastq read name pattern
    if(util::starts_with(r->name, "@NS") || 
       util::starts_with(r->name, "@NB") || 
       util::starts_with(r->name, "@A0")){
       delete r;
       return true;
    }
    
    delete r;
    return false;
}

void Evaluator::evaluateSeqLen(int& r1Len, int& r2Len){
    if(!this->r1File.empty()){
        r1Len = this->computeSeqLen(this->r1File);
    }
    if(!this->r2File.empty()){
        r2Len = this->computeSeqLen(this->r2File);
    }
}

int Evaluator::computeSeqLen(const std::string& filename){
    FqReader fqr(filename);
    size_t records = 0;
    Read* r = NULL;

    int seqLen = 0;
    while(records < 1000){
        r = fqr.read();
        if(!r){
            break;
        }
        seqLen = std::max(seqLen, r->length());
        ++records;
        delete r;
    }
    return seqLen;
}

void Evaluator::computeOverRepSeq(const std::string& filename, std::map<std::string, size_t>& hotSeqs, const int& seqLen){
    FqReader fqr(filename);
    std::map<std::string, size_t> seqCounts;
    const size_t BASE_LIMIT = 151 * 10000;
    size_t records = 0;
    size_t bases = 0;
    Read* r = NULL;
    int rlen = 0;
    std::set<int> steps = {10, 20, 40, 100, std::min(150, seqLen - 2)};
    std::string seq;
    size_t count;

    while(bases < BASE_LIMIT){
        r = fqr.read();
        if(!r){
            break;
        }
        rlen = r->length();
        bases += rlen;
        ++records;
        for(auto& s: steps){
            int step = steps[s];
            for(int i = 0; i < rlen - step; ++i){
                seq = r->seq.seqStr.substr(i, step);
                if(seqCounts.count(seq) > 0){
                    ++seqCounts[seq];
                }else{
                    seqCounts[seq] = 1;
                }
            }
        }
        delete r;
    }

    for(auto& e: seqCounts){
        seq = e.first;
        count = e.second;

        if((seq.length() >= seqLen - 1 && count >= 3) ||
           (seq.length() >= 100 && count >= 5)        ||
           (seq.length() >= 40 && count >= 20)        ||
           (seq.length() >= 20 && count >= 100)       ||
           (seq.length() >= 10 && count >= 500)){
            hotSeqs[seq] = count;
        }
    }
    
    bool isSubstring = false;
    std::string seq2;
    size_t count2 = 0;
    std::map<std::string, size_t>::iterator iter, iter2;
    iter = hotSeqs.begin();
    while(iter != hotSeqs.end()){
        seq = iter->first;
        count = iter->second;
        isSubstring = false;
        for(iter2 = hotSeqs.begin(); iter2 != hotSeqs.end(); ++iter2){
            seq2 = iter2->first;
            count2 = iter2->second;
            if(seq != seq2 && seq2.find(seq) != std::string::npos && count / count2 < 10){
                isSubstring = true;
                break;
            }
        }
        if(isSubstring){
            hotSeqs.erase(iter++);
        }else{
            ++iter;
        }
    }
}

void Evaluator::evaluateReadNum(size_t& readNum){
    FqReader fqr(this->r1File);
    const size_t READ_LIMIT = 512 * 1024;
    const size_t BASE_LIMIT = 151 * 512 * 1024;
    size_t records = 0;
    size_t bases = 0;
    size_t firstReadPos = 0;
    size_t bytesRead;
    size_t bytesTotal;
    bool reachedEOF = false;
    bool first = true;
    Read* r = NULL;
    while(records < READ_LIMIT && bases < BASE_LIMIT){
        r = fqr.read();
        if(!r){
            reachedEOF = true;
            break;
        }
        if(first){
            fqr.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        ++records;
        bases += r->length();
        delete r;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    }else if(records > 1){
        fqr.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) (records - 1);
        // bytesPerRead is a little over estimated due to some reasons unknown
        readNum = (size_t)(bytesTotal * 1.01 / bytesPerRead);
    }
}

std::string Evaluator::evalAdapterAndReadNum(size_t& readNum, bool isR2, size_t trim){
    std::string filename = (isR2 ? this->r2File : this->r1File);
    FqReader fqr(filename);
    const size_t READ_LIMIT = 256 * 1024;
    const size_t BASE_LIMIT = 151 * READ_LIMIT;
    size_t records = 0;
    size_t bases = 0;
    size_t firstReadPos = 0;
    size_t bytesRead;
    size_t bytesTotal;

    Read** loadedReads = new Read*[READ_LIMIT];
    std::memset(loadedReads, 0, sizeof(Read*) * READ_LIMIT);
    bool reachedEOF = false;
    bool first = true;
    size_t rlen = 0;

    while(records < READ_LIMIT && bases < BASE_LIMIT){
        Read* r = fqr.read();
        if(!r){
            reachedEOF = true;
            break;
        }
        if(first){
            fqr.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        bases += r->length();
        loadedReads[records] = r;
        ++records;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    }else if(records > 1){
        fqr.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (records - 1);
        readNum = (size_t)(bytesTotal * 1.01 / bytesPerRead);
    }

    if(records < 10000){
        for(size_t i = 0; i < records; ++i){
            delete loadedReads[i];
            loadedReads[i] = NULL;
        }
        return "";
    }
    
    // shit the last cycle(s) for evaluation since it is so noisy, especially for illumina data
    const size_t shiftTail = std::max(1, trim);
    const size_t keylen = 10;
    size_t size = 1 << (keylen * 2);
    size_t* counts = new size_t[size];
    std::memset(counts, 0, sizeof(size_t) * size);
    Read* r = NULL;
    int k = -1;
    for(size_t i = 0; i < records; ++i){
        r = loadedReads[i];
        const char* data = r->seq.seqStr.c_str();
        key = -1;
        for(size_t pos = 20; pos <= r->length() - keylen - shiftTail; ++pos){
            key = this->seq2int(r->seq.seqStr, pos, keylen, key);
            if(key >= 0){
                ++counts[key];
            }
        }
    }
    // set AAAAAAAAAA # = 0
    counts[0] = 0;
    // get the top N
    const int topnum = 10;
    int topkeys[topnum] = {0};
    size_t total = 0;
    int atcg[4] = {0};
    int baseOfBit = 0;
    bool lowComplexity = false;
    for(size_t k = 0; k < size; ++k){
        atcg[4] = {0};
        for(size_t i = 0; i < keylen; ++i){
            baseOfBit = (k >> (i * 2)) & 0x03;
            ++atcg[baseOfBit];
        }
        lowComplexity = false;
        // keylen = 10, if one kind of base have at least 6 in this sequence, then it is a lowComplexity seqence
        for(int b = 0; b < 4; ++b){
            if(atcg[b] >= keylen - 4){
                lowComplexity = true;
                break;
            }
        }
        if(lowComplexity){
            continue;
        }
        // keylen = 10, if GC counts is at least 8 in this sequence, then it is a high gc content sequence
        if(atcg[2] + atcg[3] >= keylen - 2){
            continue;
        }
        // if sequence starts with GGGG, also skip
        if(k >> 12 == 0xff){
            continue;
        }

        size_t val = counts[k];
        total += val;
        for(int t = topnum - 1; t >= 0; --t){
            if(val < counts[topkeys[t]]){
                // shit topkeys[t + 2] ... topkeys[topnum - 1] and insert
                if(t < topnum - 1){
                    for(int m = topnum - 1; m > t + 1; --m){
                        topkeys[m] = topkeys[m - 1];
                    }
                    topkeys[t+1] = k;
                }
                break;
            }else if(t == 0){
                // find the max ever since
                for(int m = topnum - 1; m > t; --m){
                    topkeys[m] = topkeys[m - 1];
                }
                topkeys[t] = k;
            }
        }
    }

    const int FOLD_THRESHOLD = 20;
    for(int t = 0; t < topnum; ++t){
        int key = topkeys[t];
        if(key == 0){
            continue;
        }
        std::string seq = this->int2seq(key, keylen);
        size_t count = counts[key];
        // if count < 10 || count < n * p ( n = FOLD_THRESHOLD, p = total/size)
        if(count < 10 || count * size < total * FOLD_THRESHOLD){
            break;
        }
        // skip low complexity seq, this may not be needed as low complexity tested before contains this
        int diff = 0;
        for(int s = 0; s < seq.length - 1; ++s){
            if(seq[s] != seq[s+1]){
                ++diff;
            }
        }
        if(diff < 3){
            continue;
        }
        std::string adapter = this->getAdapterWithSeed(key, loadedReads, records, keylen);
        if(!adapter.empty()){
            delete[] counts;
            for(int r = 0; r < records; ++r){
                delete loadedReads[r];
                loadedReads[r] = NULL;
            }
            delete[] loadedReads;
            return adapter;
        }
    }
            
    delete[] counts;
    for(int r = 0; r < records; ++r){
        delete loadedReads[r];
        loadedReads[r] = NULL;
    }
    delete[] loadedReads;
    return "";
}

std::string Evaluator::getAdapterWithSeed(int seed, Read** loadedReads, long records, int keylen, int trim){
    const int shiftTail = std::max(1, trim);
    NucleotideTree forwardTree, backwardTree;
    //construct trees
    for(int i = 0; i < records; ++i){
        Read* r = loadedReads[i];
        const char* data = r->seq.seqStr.c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length() - keylen - shiftTail; ++pos){
            key = seq2int(r->seq.seqStr, pos, keylen, key);
            if(key == seed){
                forwardTree.addSeq(r->seq.seqStr.substr(pos + keylen, r->length() - keylen - shiftTail - pos));
                backwardTree.addSeq(util::reverse(r->seq.seqStr.substr(0, pos)));
            }
        }
    }
    bool reachedLeaf = true;
    std::string forwardPath = forwardTree.getDominantPath(reachedLeaf);
    std::string backwardPath = backwardTree.getDominantPath(reachedLeaf);
    std::string adapter = util::reverse(backwardPath) + this->int2seq(seed, keylen) + forwardPath;
    if(adapter.length() > 60){
        adapter.resize(60);
    }
    std::string matchedAdapter = this->matchKnownAdapter(adapter);
    if(!matchedAdapter.empty()){
        std::map<std::string, std::string> knownAdapters = ::getKnownAdapter();
        std::cerr << knownAdapters[matchedAdapter] << ": " << matchedAdapter << std::endl;
        return matchedAdapter;
    }else{
        if(reachedLeaf){
            std::cerr << adapter << endl;
            return adapter;
        }else{
            return "";
        }
    }
}

std::string Evaluator::matchKnownAdapter(std::string seq){
    std::map<std::string, std::string> knownAdapters = ::getKnownAdapter();
    std::map<std::string, std::string>::iterator iter;
    for(iter = knownAdapters.begin(); iter != knownAdapters.end(); ++iter) {
        std::string adapter = iter->first;
        std::string desc = iter->second;
        if(seq.length() < adapter.length()) {
            continue;
        }
        int diff = 0;
        for(int i=0; i < adapter.length() && i < seq.length(); ++i) {
            if(adapter[i] != seq[i])
                ++diff;
        }
        if(diff == 0)
            return adapter;
    }
    return "";
}

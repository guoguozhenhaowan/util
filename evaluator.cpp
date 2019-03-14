#include "evaluator.h"
#include "knownadapters.h"

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
}





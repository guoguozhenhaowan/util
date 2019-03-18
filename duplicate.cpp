#include "duplicate.h"

namespace fqlib{
    Duplicate::Duplicate(Options* opt){
        mOptions = opt;
        mKeyLenInBase = opt->duplicate.keylen;
        mKeyLenInBit = 1 << (2 * mKeyLenInBase);
        mDups = new u_int64_t[mKeyLenInBit];
        std::memset(mDups, 0, sizeof(u_int64_t) * mKeyLenInBit);
        mCounts = new u_int32_t[mKeyLenInBit];
        std::memset(mCounts, 0, sizeof(u_int32_t) * mKeyLenInBit);
        mGC = new u_int8_t[mKeyLenInBit];
        std::memset(mGC, 0, sizeof(u_int8_t) * mKeyLenInBit);
    }

    Duplicate::~Duplicate(){
        delete[] mDups;
        delete[] mCounts;
        delete[] mGC;
    }

    u_int64_t Duplicate::seq2int(const char* cstr, int start, int keylen, bool& valid){
        u_int64_t ret = 0;
        for(int i = 0; i < keylen; ++i){
            ret <<= 2;
            switch(cstr[start + i]){
                case 'A':
                    ret += 0;
                    break;
                case 'T':
                    ret += 1;
                    break;
                case 'C':
                    ret += 2;
                    break;
                case 'G':
                    ret += 3;
                    break;
                default:
                    valid = false;
                    return 0;
            }
        }
        return ret;
    }

    void Duplicate::addRecord(u_int32_t key, u_int64_t kmer32, u_int8_t gc){
        if(mCounts[key] == 0){
            mCounts[key] = 1;
            mDups[key] = kmer32;
            mGC[key] = gc;
        }else{
            if(mDups[key] == kmer32){
                ++mCounts[key];
            }else if(mDups[key] > kmer32){
                mDups[key] = kmer32;
                mCounts[key] = 1;
                mGC[key] = gc;
            }
        }
    }

    void Duplicate::statRead(Read* r){
        if(r->length() < 32){
            return;
        }

        int start1 = 0;
        int start2 = std::max(0, r->length() - 32 - 5);

        const char* cstr = r->seq.seqStr.c_str();
        bool valid = true;
        u_int64_t ret = seq2int(cstr, start1, mKeyLenInBase, valid);
        u_int32_t key = (u_int32_t)ret;
        if(!valid){
            return;
        }
        u_int64_t kmer32 = seq2int(cstr, start2, 32, valid);
        if(!valid){
            return;
        }
        u_int8_t gc = 0;

        if(mCounts[key] == 0){
            for(int i = 0; i < r->length(); ++i){
                if(cstr[i] == 'C' || cstr[i] == 'G'){
                    ++gc;
                }
            }
        }
        gc = std::round(255.0 * (double) gc / (double) r->length());
        addRecord(key, kmer32, gc);
    }

    void Duplicate::statPair(Read* r1, Read* r2){
        if(r1->length() < 32 || r2->length() < 32){
            return;
        }
        const char* cstr1 = r1->seq.seqStr.c_str();
        const char* cstr2 = r2->seq.seqStr.c_str();
        bool valid = true;

        u_int64_t ret = seq2int(cstr1, 0, mKeyLenInBase, valid);
        u_int32_t key = (u_int32_t)ret;
        if(!valid){
            return;
        }
        u_int64_t kmer32 = seq2int(cstr2, 0, 32, valid);
        if(!valid){
            return;
        }

        u_int8_t gc = 0;
        if(mCounts[key] == 0){
            for(int i = 0; i < r1->length(); ++i){
                if(cstr1[i] == 'C' || cstr1[i] == 'G'){
                    ++gc;
                }
            }
            for(int i = 0; i < r2->length(); ++i){
                if(cstr2[i] == 'C' || cstr2[i] == 'G'){
                    ++gc;
                }
            }
        }

        gc = std::round(255.0 * (double) gc / (double)(r1->length() + r2->length()));

        addRecord(key, kmer32, gc);
    }

    double Duplicate::statAll(int* hist, double* meanGC, int histSize){
        size_t totalNum = 0;
        size_t dupNum = 0;
        int* gcStatNum = new int[histSize];
        std::memset(gcStatNum, 0, sizeof(int) * histSize);
        for(int key = 0; key < mKeyLenInBit; ++key){
            u_int32_t count = mCounts[key];
            u_int8_t gc = mGC[key];
            if(count > 0){
                totalNum += count;
                dupNum += count - 1;
                if(count > histSize){
                    ++hist[histSize -1];
                    meanGC[histSize -1] += gc;
                    ++gcStatNum[histSize -1];
                }else{
                    ++hist[count];
                    meanGC[histSize -1] += gc;
                    ++gcStatNum[histSize -1];
                }
            }
        }
        
        for(int i = 0; i < histSize; ++i){
            if(gcStatNum[i] > 0){
                meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
            }
        }

        delete[] gcStatNum;

        if(totalNum == 0){
            return 0.0;
        }else{
            return (double)dupNum / (double)totalNum;
        }
    }
}

#include "polyx.h"

namespace fqlib{
    PolyX::PolyX(){
    }

    PolyX::~PolyX(){
    }

    void PolyX::trimPolyG(Read* r1, Read* r2, FilterResult* fr, int compareReq){
        trimPolyG(r1, fr, compareReq);
        trimPolyG(r2, fr, compareReq);
    }

    void PolyX::trimPolyG(Read* r, FilterResult* fr, int compareReq){
        const int allowedOneMismatchForEach = 8;
        const int maxMismatch = 5;
        const char* data = r->seq.seqStr.c_str();
        int rlen = r->length();
        int mismatch = 0;
        int i = 0;
        int firstGpos = rlen - 1;
        for(int i = 0; i < rlen; ++i){
            if(data[rlen - i - 1] != 'G'){
                ++mismatch;
            }else{
                firstGpos = rlen - i - 1;
            }
            int allowedMismatch = (i + 1)/allowedOneMismatchForEach;
            if(mismatch > maxMismatch || (mismatch > allowedMismatch && i >= compareReq - 1)){
                break;
            }
        }
        if(i >= compareReq){
            r->resize(firstGpos);
        }
    }

    void PolyX::trimPolyX(Read* r1, Read* r2, FilterResult* fr, int compareReq){
        trimPolyX(r1, fr, compareReq);
        trimPolyX(r2, fr, compareReq);
    }

    void PolyX::trimPolyX(Read* r, FilterResult* fr, int compareReq){
        const int allowedOneMismatchForEach = 8;
        const int maxMismatch = 5;
        const char* data = r->seq.seqStr.c_str();
        int rlen = r->length();
        int atcgNumbers[4] = {0, 0, 0, 0};
        char atcgBases[4] = {'A', 'T', 'C', 'G'};
        int pos = 0;
        for(pos = 0; pos < rlen; ++pos){
            switch(data[rlen - 1 -pos]){
                case 'A':
                    ++atcgNumbers[0];
                    break;
                case 'T':
                    ++atcgNumbers[1];
                    break;
                case 'C':
                    ++atcgNumbers[2];
                    break;
                case 'G':
                    ++atcgNumbers[3];
                    break;
                case 'N':
                    for(int i = 0; i < 4; ++i){
                        ++atcgNumbers[i];
                    }
                    break;
                default:
                    break;
            }
            int cmp = (pos + 1);
            int allowedMismatch = std::min(maxMismatch, cmp/allowedOneMismatchForEach);
            bool needToBreak = true;
            for(int b = 0; b < 4; ++b){
                if(cmp - atcgNumbers[b] <= allowedMismatch){
                    needToBreak = false;
                }
            }
            if(needToBreak && (pos >= allowedOneMismatchForEach || pos + 1 >= compareReq - 1)){
                break;
            }
        }

        if(pos + 1 >- compareReq){
            int poly;
            int maxCount = -1;
            for(int b = 0; b < 4; ++b){
                if(atcgNumbers[b] > maxCount){
                    maxCount = atcgNumbers[b];
                    poly = b;
                }
            }
            char polyBase = atcgBases[poly];
            while(data[rlen - pos - 1] != polyBase && pos >= 0){
                --pos;
            }
            r->resize(rlen - pos - 1);
        }
    }
}

#include "filter.h"

int Filter::passFilter(Read* r){
    if(r == NULL || r->length() == 0){
        return compar::FAIL_LENGTH;
    }

    int rlen = r->length();
    int lowQualNum = 0;
    int nBaseNum = 0;
    const char* seq = r->seq.seqStr.c_str();
    const char* qual = r->quality.c_str();

    if(this->lowQualThreshold > 0 || this->maxNBaseNum > 0 || this->minLen > 0 || this->maxLen > 0){
        for(int i = 0; i < rlen; ++i){
            if(seq[i] == 'N'){
                ++nBaseNum;
            }
            if(qual[i] < this->lowQualThreshold){
                ++lowQualNum;
            }
        }
    }

    if(this->lowQualThreshold > 0 && lowQualNum > this->maxLowQualNum){
        return compar::FAIL_QUALITY;
    }

    if(this->maxNBaseNum > 0 && nBaseNum > this->maxNBaseNum){
        return compar::FAIL_N_BASE;
    }

    if(this->minLen > 0 && rlen < this->minLen){
        return compar::FAIL_LENGTH;
    }

    if(this->maxLen > 0 && rlen > this->maxLen){
        return compar::FAIL_TOO_LONG;
    }

    if(this->minComplexity > 0 && !(this->passLowComplexityFliter(r))){
        return compar::FAIL_COMPLEXITY;
    }

    return compar::PASS_FILTER;
}

bool Filter::passLowComplexityFliter(Read* r){
    int diff = 0;
    int rlen = r->length();
    if(rlen <= 1){
        return false;
    }
    const char* seq = r->seq.seqStr.c_str();
    for(int i = 0; i < rlen - 1; ++i){
        if(seq[i] != seq[i+1]){
            ++diff;
        }
    }
    return (double)diff/(rlen - 1) >= this->minComplexity;
}

Read* Filter::trimAndCut(Read* r, int forceFrontCut, int forceTailCut, QualCutOpt* opt){
    // do not need quality cutting
    if(forceFrontCut == 0 && forceTailCut == 0 && !opt->cutFront && !opt->cutTail && !opt->cutRright){
        return r;
    }

    int rlen = r->length() - forceFrontCut - forceTailCut;
    if(rlen < 0){
        return NULL;
    }

    if(forceFrontCut == 0 && !opt->cutFront && !opt->cutTail && !opt->cutRright){
        r->resize(rlen);
        return r;
    }else if(!opt->cutFront && !opt->cutTail && !opt->cutRright){
        r->seq.seqStr = r->seq.seqStr.substr(forceFrontCut, rlen);
        r->quality = r->quality.substr(forceFrontCut, rlen);
        return r;
    }

    // need quality cutting
    int l = r->length();
    const char* qual = r->quality.c_str();
    const char* seq = r->seq.seqStr.c_str();
    // quality cutting forward by sliding window
    if(opt->cutFront){
        int w = opt->windowFront;
        int s = forceFrontCut;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[s + i];
        }
        for(s = forceFrontCut; s + w < l - forceTailCut; ++s){
            totalQual += qual[s + w -1];
            if(s > forceFrontCut){
                totalQual -= qual[s - 1];
            }
            if((double)totalQual/(double)w >= 33 + opt->minFrontQual){
                break;
            }
        }
        if(s > 0){
            s = s + w -1;
        }
        // forceFrontCut is relocated one base before the last sliding window's last base and skip any N afterwards aswell
        while(s < l && seq[s] == 'N'){
            ++s;
        }
        forceFrontCut = s;
        rlen = l - forceFrontCut - forceTailCut;
    }

    // quality cutting from right by sliding window
    if(opt->cutRright){
        int w = opt->windowRight;
        int s = forceFrontCut;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        bool foundLowQualWindow = false;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[s + i];
        }
        for(s = forceFrontCut; s + w < l - forceTailCut; ++s){
            totalQual += qual[s + w -1];
            if(s > forceFrontCut){
                totalQual -= qual[s - 1];
            }
            if((double)totalQual/(double)w < 33 + opt->minRightQual){
                foundLowQualWindow = true;
                break;
            }
        }
        if(foundLowQualWindow){
            while(s < l - 1 && qual[s] >= 33 + opt->minRightQual){
                ++s;
            }
            rlen = s - forceFrontCut;
        }
    }

    // quality cutting backward by sliding window
    if(!opt->cutRright && opt->cutTail){
        int w = opt->windowTail;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        int t = l - forceTailCut - 1;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[t - i];
        }
        for(t = l - forceTailCut - 1; t - w >= forceFrontCut; --t){
            totalQual += qual[t - w + 1];
            if(t < l - forceTailCut - 1){
                totalQual -= qual[t + 1];
            }
            if((double)totalQual/(double)w >= 33 + opt->minTailQual){
                break;
            }
        }
        if(t < l - 1){
            t = t - w + 1;
        }
        while(t >= 0 && seq[t] == 'N'){
            --t;
        }
        rlen = t - forceFrontCut + 1;
    }

    if(rlen <= 0 || forceFrontCut >= l - 1){
        return NULL;
    }
    r->seq.seqStr = r->seq.seqStr.substr(forceFrontCut, rlen);
    r->quality = r->quality.substr(forceFrontCut, rlen);
    return r;
}

bool Filter::match(const std::vector<std::string>& list, const std::string& target, int threshold){
    int tlen = target.length();
    int slen = 0;
    int diff = 0;
    for(int i = 0; i < list.size(); ++i){
        diff = 0;
        slen = list[i].length();
        for(int s = 0; s < slen && s < tlen; ++s){
            if(list[i][s] != target[s]){
                ++diff;
                if(diff > threshold){
                    break;
                }
            }
        }
        if(diff <= threshold){
            return true;
        }
    }
    return false;
}

#include "stats.h"

Stats::Stats(const std::string& fileName){
    this->fqFile = fileName;
    Evaluator e(fileName, 0);
    int readLen = e.getReadLen();
    this->reads = 0;
    this->bases = 0;
    this->evaluatedSeqLen = readLen;
    this->cycles = readLen;
    this->bufLen = readLen + 1024;
    this->q20Total = 0;
    this->q30Total = 0;
    this->summarized = false;
    this->kmerMax = 0;
    this->kmerMin = 0;
    this->kmerLen = 5;
    this->kmerBufLen = 2 << (this->kmerLen * 2);
    this->lengthSum = 0;
    this->overRepSampleFreq = 100;
}

void Stats::allocateRes(){
    for(int i = 0; i < 8; ++i){
        this->q20Bases[i] = 0;
        this->q30Bases[i] = 0;
        this->baseContents[i] = 0;
        this->cycleQ20Bases[i] = new size_t[this->bufLen];
        std::memset(this->cycleQ20Bases[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleQ30Bases[i] = new size_t[this->bufLen];
        std::memset(this->cycleQ30Bases[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleBaseContents[i] = new size_t[this->bufLen];
        std::memset(this->cycleBaseContents[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleBaseQual[i] = new size_t[this->bufLen];
        std::memset(this->cycleBaseQual[i], 0, sizeof(size_t) * this->bufLen);
    }

    this->cycleTotalBase = new size_t[this->bufLen];
    std::memset(this->cycleTotalBase, 0, sizeof(size_t) * this->bufLen);
    this->cycleTotalQual = new size_t[this->bufLen];
    std::memset(this->cycleTotalQual, 0, sizeof(size_t) * this->bufLen);
    
    if(this->kmerLen){
        this->kmer = new size_t[this->kmerBufLen];
        std::memset(this->kmer, 0, sizeof(size_t) * this->kmerBufLen);
    }
    
    if(this->overRepSampleFreq){
        this->initOverRepSeq();
    }
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= this->bufLen){
        return;
    }
    size_t* newBuf = NULL;

    for(int i = 0; i < 8; ++i){
        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleQ20Bases[i], sizeof(size_t) * this->bufLen);
        delete this->cycleQ20Bases[i];
        this->cycleQ20Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleQ30Bases[i], sizeof(size_t) * this->bufLen);
        delete this->cycleQ30Bases[i];
        this->cycleQ30Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleBaseContents[i], sizeof(size_t) * this->bufLen);
        delete this->cycleBaseContents[i];
        this->cycleBaseContents[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleBaseQual[i], sizeof(size_t) * this->bufLen);
        delete this->cycleBaseQual[i];
        this->cycleBaseQual[i] = newBuf;
    }

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, this->cycleTotalBase, sizeof(size_t) * this->bufLen);
    delete this->cycleTotalBase;
    this->cycleTotalBase = newBuf;

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, this->cycleTotalQual, sizeof(size_t) * this->bufLen);
    delete this->cycleTotalQual;
    this->cycleTotalQual = newBuf;

    this->bufLen = newBufLen;
}

Stats::~Stats(){
    for(int i = 0; i < 8; ++i){
        delete cycleQ20Bases[i];
        cycleQ20Bases[i] = NULL;
        
        delete cycleQ30Bases[i];
        cycleQ30Bases[i] = NULL;
        
        delete cycleBaseContents[i];
        cycleBaseContents[i] = NULL;
        
        delete cycleBaseQual[i];
        cycleBaseQual[i] = NULL;
    }
    
    delete cycleTotalBase;
    delete cycleTotalQual;
    
    for(auto& e : this->qualityCurves){
        delete e.second;
    }
    
    for(auto& e : this->contentCurves){
        delete e.second;
    }
    
    if(this->kmerLen){
        delete this->kmer;
    }

    deleteOverRepSeqDist();
}

void Stats::summarize(bool forced){
    if(this->summarized && !forced){
        return;
    }

    // cycle and total bases
    int c = 0;
    for(c = 0; c < this->bufLen; ++c){
        this->bases += this->cycleTotalBase[c];
        if(this->cycleTotalBase[c] == 0){
            break;
        }
    }
    this->cycles = c;

    // Q20, Q30, base content
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < this->cycles; ++j){
            this->q20Bases[i] += this->cycleQ20Bases[i][j];
            this->q30Bases[i] += this->cycleQ30Bases[i][j];
            this->baseContents[i] += this->cycleBaseContents[i][j];
        }
        this->q20Total += this->q20Bases[i];
        this->q30Total += this->q30Bases[i];
    }

    // quality curve for mean qual
    double* meanQualCurve = new double[this->cycles];
    std::memset(meanQualCurve, 0, sizeof(double) * this->cycles);
    for(int i = 0; i < this->cycles; ++i){
        meanQualCurve[i] = (double)this->cycleTotalQual[i] / (double)this->cycleTotalBase[i];
    }
    this->qualityCurves["mean"] = meanQualCurve;

    // quality curves and base contents curves for different nucleotides
    char nucleotides[5] = {'A', 'T', 'C', 'G', 'N'};
    for(int i = 0; i < 5; ++i){
        char b = nucleotides[i] & 0x07;
        double* qualCurve = new double[this->cycles];
        std::memset(qualCurve, 0, sizeof(double) * this->cycles);
        double* contentCurve = new double[this->cycles];
        std::memset(contentCurve, 0, sizeof(double) * this->cycles);
        for(int j = 0; j < this->cycles; ++j){
            qualCurve[j] = (double)this->cycleBaseQual[b][j] / (double)this->cycleBaseContents[b][j];
            contentCurve[j] = (double)this->cycleBaseContents[b][j] / (double)this->cycleBaseContents[b][j];
        }
        this->qualityCurves[std::string(1, nucleotides[i])] = qualCurve;
        this->contentCurves[std::string(1, nucleotides[i])] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[this->cycles];
    std::memset(gcContentCurve, 0, sizeof(double) * this->cycles);
    char gIndex = 'G' & 0x07;
    char cIndex = 'C' & 0x07;
    for(int i = 0; i < this->cycles; ++i){
        gcContentCurve[i] = (double)(this->cycleBaseContents[gIndex][i] + this->cycleBaseContents[cIndex][i])/ (double)this->cycleTotalBase[i];
    }
    this->contentCurves["GC"] = gcContentCurve;

    if(this->kmerLen){
        // K-mer statistics
        this->kmerMin = this->kmer[0];
        this->kmerMax = this->kmer[0];
        for(int i = 1; i < this->kmerBufLen; ++i){
            this->kmerMin = std::min(this->kmer[i], this->kmerMin);
            this->kmerMax = std::max(this->kmer[i], this->kmerMax);
        }
    }

    this->summarized = true;
}

int Stats::getMeanLength(){
    if(this->reads = 0){
        return 0;
    }
    return this->lengthSum/this->reads;
}

void Stats::statRead(Read* r){
    int len = r->length();
    this->lengthSum += len;
    if(this->bufLen < len){
        this->extendBuffer(std::max(len + 100, (int)(len * 1.5)));
    }
    const char* seqStr = r->seq.seqStr.c_str();
    const char* qualStr = r->quality.c_str();
    int kmerVal = 0;
    bool needFullCompute = true; //recompute kmerVal needed
    int kmerMask = ~(0x3);
    for(int i = 0; i < len; ++i){
        char base = seqStr[i];
        char qual = qualStr[i];
        char bIndex = base & 0x07;
        const char q20 = '5';
        const char q30 = '?';
        if(qual > q20){
            ++(this->cycleQ20Bases[bIndex][i]);
        }
        if(qual > q30){
            ++this->cycleQ20Bases[bIndex][i];
            ++this->cycleQ30Bases[bIndex][i];
        }
        ++this->cycleBaseContents[bIndex][i];
        this->cycleBaseQual[bIndex][i] += (qual - 33);
        ++this->cycleTotalBase[i];
        this->cycleTotalQual[i] += (qual - 33);
        if(base == 'N'){
            needFullCompute = true;
            continue;
        }

        if(this->kmerLen){
            if(i < this->kmerLen - 1){
                continue;
            }
            if(!needFullCompute){
                int val = this->base2val(base);
                if(val < 0){
                    needFullCompute = true;
                    continue;
                }else{
                    kmerVal = ((kmerVal << 2) & kmerMask) | val;
                    ++this->kmer[kmerVal];
                }
            }else{
                bool valid = true;
                kmerVal = 0;
                for(int k = 0; k < this->kmerLen; ++k){
                    int val = this->base2val(this->seqStr[i - this->kmerLen + k]);
                    if(val < 0){
                        valid = false;
                        break;
                    }
                    kmerVal = ((kmerVal << 2) & kmerMask) | val;
                }
                if(!valid){
                    needFullCompute = true;
                    continue;
                }else{
                    ++this->kmer[kmerVal];
                    needFullCompute = false;
                }
            }
        }

        if(this->overRepSampleFreq){
            if(this->reads % this->overRepSampleFreq == 0){
                std::set<int> steps[5] = {10, 20, 40, 100, std::min(150, this->evaluatedSeqLen - 2)};
                for(auto s: steps){
                    for(int i=0; i< len - step; ++i){
                        std::string seq = r->seq.seqStr.substr(i, step);
                        if(this->overRepSeq.count(seq) > 0){
                            ++this->overRepSeq.count(seq);
                        }else{
                            for(int p = i; p < seq.length() + i && p < this->evaluatedSeqLen; ++p){
                                ++this->overRepSeqDist[seq][p];
                            }
                        }
                        i += step;
                    }
                }
            }
        }
    }
    ++this->reads;
}

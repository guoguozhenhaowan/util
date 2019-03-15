#ifndef STATS_H
#define STATS_H

#include <map>
#include <memory>
#include <vector>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <sstream>
#include "read.h"
#include "util.h"
#include "evaluator.h"

/** Class to do statistics of a fastq file */
class Stats{
    public:
        std::string fqFile;            ///< fastq file name to be analysised
        size_t reads;                  ///< total reads of this->fqFile
        size_t bases;                  ///< total bases of this->fqFile
        int evaluatedSeqLen;           ///< estimated read length of this->fqFile
        int cycles;                    ///< maximum cycle of this->fqFile
        int bufLen;                    ///< buffer length to store a static item of one read in this->fqFile
        size_t q20Total;               ///< total bases with quality equal or greater than 20 in this->fqFile
        size_t q30Total;               ///< total bases with quality equal or greater than 30 in this->fqFile
        bool summarized;               ///< summarized or not, if this->summarize() called, then this->summarized = true
        size_t kmerMax;                ///< the maximum kmerLen kmer in this->fqFile
        size_t kmerMin;                ///< the minimum kmerLen kmer in this->fqFile
        size_t kmerLen = 5;            ///< kmer length to calculate, default 5
        int kmerBufLen;                ///< kmer statistic array length, just equal 2 << (2 * this->kmerLen)
        size_t lengthSum;              ///< total length of reads in this->fqFile
        size_t q20Bases[8];            ///< each base(ATCG) counts with quality equal or greater than 20 in this->fqFile
        size_t q30Bases[8];            ///< each base(ATCG) counts with quality equal or greater than 30 in this->fqFile
        size_t baseContents[8];        ///< each base(ATCG) counts of this->fqFile
        size_t *cycleQ20Bases[8];      ///< each base(ATCG) counts at each cycle with quality equal or greater than 20 in this->fqFile
        size_t *cycleQ30Bases[8];      ///< each base(ATCG) counts at each cycle with quality equal or greater than 30 in this->fqFile
        size_t *cycleBaseContents[8];  ///< each base(ATCG) counts at each cycle of this->fqFile
        size_t *cycleBaseQual[8];      ///< each base(ATCG) quality at each cycle of this->fqFile
        size_t *cycleTotalBase;        ///< total base counts of each cycle
        size_t *cycleTotalQual;        ///< total quality of each cycle
        size_t *kmer;                  ///< kmer count array 
        int overRepSampleFreq = 100;   ///< over representation analysis sampling frequence, default 100
         

        std::map<std::string, double*> qualityCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle quality 
        std::map<std::string, double*> contentCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle content
        std::map<std::string, size_t> overRepSeq;
        std::map<std::string, size_t*> overRepSeqDist;

    public:
        /** Construct a Stats object of fastq fileName
         * Get estimated read length
         */
        Stats(const std::string& fileName);
        
        /** Destroy a Stats object
         * Free allocated mamory
         */
        ~Stats();

        /** Set Kmer length to calculate
         * @param kLen Kmer length, if kLen = 0, Kmer statistics will not be done
         */ 
        inline void setKmerLen(const int kLen){
            this->kmerLen = kLen;
        }

        /** Set over representation analysis sampling frequence
         * @param readFreq if readFreq = 0, over representation analysis will not be done
         */
        inline void setOverRepSampleFreq(int readFreq){
            this->overRepSampleFreq = readFreq;
        }

        /** Allocate resources for a Stats object
         */
        void allocateRes(){
        }

        int getCycles();
        size_t getReads();
        size_t getBases();
        size_t getQ20();
        size_t getQ30();
        size_t getGCNumber();
        void statRead(Read* r);

        static Stats* merge(std::vector<Stats*>& list);
        void print();

        /** Summary statistic items, get this->bases, this->cycles, this->q20Bases, this->q30Bases, this->q20Total, this->q30Total,
         * this->qualityCurves["mean"], this->qualityCurves["A/T/C/G/N"], this->contentCurves["A/T/C/G/N"], this->contentCurves["GC"]
         * this->kmerMin, this->kmerMax, after summary, this->summarize will be set to true
         * @param forced forced to run summarize
         * this->summarize will only run if (this->summarized is false) || (this->summarized is true && force = false)
         */
        void summarize(bool forced = false);

        void reportJson(std::ofstream& ofs, std::string padding);
        void reportHtml(std::ofstream& ofs, std::string filteringType, std::string readName);
        void reportHtmlQuality(std::ofstream& ofs, std::string filteringType, std::string readName);
        void reportHtmlKmer(std::ofstream& ofs, std::string filteringType, std::string readName);
        void reportHtmlORA(std::ofstream& ofs, std::string filteringType, std::string readName);
        bool isLongRead();
        void initOverRepSeq();
        int getMeanLength();

    public:
        static std::string list2string(double* list, int size);
        static std::string list2string(double* list, int size, size_t* coords);
        static std::string list2string(size_t* list, int size);
        static int base2val(char base);

    public:
        void extendBuffer(int newBufLen);
        std::string makeKmerTD(int i, int j);
        std::string kmer3(int val);
        std::string kmer2(int val);
        void deleteOverRepSeqDist();
        bool overRepPassed(std::string& seq, size_t count);
};

#endif

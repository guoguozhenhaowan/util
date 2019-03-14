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

class Stats{
    public:
        bool isRead2;
        size_t reads;
        int evaluatedSeqLen;
        size_t *cycleQ30Bases[8];
        size_t *cycleQ20Bases[8];
        size_t *cycleBaseContents[8];
        size_t *cycleBaseQual[8];
        size_t *cycleTotalBase;
        size_t *cycleTotalQual;
        size_t *kmer;

        std::map<std::string, double*> qualityCurves;
        std::map<std::string, double*> contentCurves;
        std::map<std::string, size_t> overRepSeq;
        std::map<std::string, size_t*> overRepSeqDist;

        int cycles;
        int bufLen;
        size_t bases;
        size_t q20Bases[8];
        size_t q30Bases[8];
        size_t baseContents[8];
        size_t q20Total;
        size_t q30Total;
        bool summarized;
        size_t kmerMax;
        size_t kmerMin;
        int kmerBufLen;
        size_t lengthSum;

    public:
        Stats(bool isR2 = false, int guessedCycles = 0, int bufferMargin = 1024);
        ~Stats();
        int getCycles();
        size_t getReads();
        size_t getBases();
        size_t getQ20();
        size_t getQ30();
        size_t getGCNumber();
        void statRead(Read* r);

        static Stats* merge(std::vector<Stats*>& list);
        void print();
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

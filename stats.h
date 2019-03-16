#ifndef STATS_H
#define STATS_H

#include <map>
#include <iomanip>
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
        size_t reads;                  ///< total reads
        size_t bases;                  ///< total bases
        int minReadLen;                ///< minimal read length
        int maxReadLen;                ///< maximum read length
        int minQual;                   ///< minimum base quality
        int maxQual;                   ///< maximum base quality
        int evaluatedSeqLen;           ///< estimated read length
        int cycles;                    ///< maximum cycle
        int bufLen;                    ///< buffer length to store a static item of one read
        size_t q20Total;               ///< number of bases with quality greater than 20 
        size_t q30Total;               ///< number of bases with quality greater than 30 
        bool summarized;               ///< summarized or not, if this->summarize() called, then this->summarized = true
        size_t kmerMax;                ///< the maximum kmerLen kmer
        size_t kmerMin;                ///< the minimum kmerLen kmer 
        int kmerLen = 5;               ///< kmer length to calculate, default 5
        int kmerBufLen;                ///< kmer statistic array length, just equal 2 << (2 * this->kmerLen)
        size_t lengthSum;              ///< total length of reads
        size_t q20Bases[8];            ///< each base(ATCG) counts with quality greater than 20
        size_t q30Bases[8];            ///< each base(ATCG) counts with quality greater than 30 
        size_t baseContents[8];        ///< each base(ATCG) counts 
        size_t *cycleQ20Bases[8];      ///< each base(ATCG) counts at each cycle with quality greater than 20 
        size_t *cycleQ30Bases[8];      ///< each base(ATCG) counts at each cycle with quality greater than 30 
        size_t *cycleBaseContents[8];  ///< each base(ATCG) counts at each cycle 
        size_t *cycleBaseQual[8];      ///< each base(ATCG) quality at each cycle 
        size_t *cycleTotalBase;        ///< total base counts of each cycle
        size_t *cycleTotalQual;        ///< total quality of each cycle
        size_t *kmer;                  ///< kmer count array 
        int overRepSampleFreq = 100;   ///< over representation analysis sampling frequence, default 100
        std::map<std::string, size_t> overRepSeq;     ///< map of <repSeq, repSeqCount>, shoulde be initialized by Evaluator before statRead
        std::map<std::string, double*> qualityCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle quality 
        std::map<std::string, double*> contentCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle content
        std::map<std::string, size_t*> overRepSeqDist;///< map of <repseq, dist*>
    
    public:
        /** Construct a Stats object of fq with length estReadLen */
        Stats(const int& estReadLen);
        
        /** Destroy a Stats object Free allocated mamory */
        ~Stats();

        /** Set Kmer length to calculate
         * @param kLen Kmer length, if kLen = 0, Kmer statistics will not be done
         */ 
        inline void setKmerLen(const int kLen){
            this->kmerLen = kLen;
            this->kmerBufLen = (1 << (this->kmerLen * 2));
        }

        /** Set over representation analysis sampling frequence
         * @param readFreq if readFreq = 0, over representation analysis will not be done
         */
        inline void setOverRepSampleFreq(int readFreq){
            this->overRepSampleFreq = readFreq;
        }

        /** Allocate resources for a Stats object */
        void allocateRes();
        
        /** Initialize overrepresentation analysis
         * @param fqEva Evaluator object  used to compute representatin sequences
         */
        void initOverRepSeq(Evaluator& fqEva);

        /** get minimum read length
         * @return minimum read length
         */
        int getMinReadLength();

        /** get maximum read length
         * @return maximum read length
         */
        int getMaxReadLength();

        /** get minimum base quality
         * @return minimum base quality
         */
        int getMinBaseQual();

        /** get maximum base quality
         * @return maximum base quality
         */
        int getMaxBaseQual();
         
        /** get cycle number
         * @return cycle number 
         */
        int getCycles();
        
        /** get total raeds number 
         * @return total read number 
         */
        size_t getReads();
        
        /** get total base number 
         * @return total base number 
         */
        size_t getBases();
        
        /** get number of bases with quality greater than 20 
         * @return number of bases with quality greater than 20
         */
        size_t getQ20();

        /** get number of bases with quality greater than 30  
         * @return number of bases with quality greater than 30 
         */
        size_t getQ30();

        /** get number of base G and C 
         * @return number of base G and C 
         */
        size_t getGCNumber();

        /** do statistics of one read 
         * @param r pointer to object Read
         */
        void statRead(Read* r);

        /** merge a list of Stats objects into one
         * @param list a list of Stats objects
         * @return a merged Stats object
         */
        static Stats* merge(std::vector<Stats*>& list);
        
        /** output Stats object to std::ostream
         * @param os std::ostream object
         * @param s Stats object
         * @return std::ostream object
         */
        friend std::ostream& operator<<(std::ostream& os, const Stats& s);

        /** Summary statistic items, get this->bases, this->cycles, this->q20Bases, this->q30Bases, this->q20Total, this->q30Total,
         * this->qualityCurves["mean"], this->qualityCurves["A/T/C/G/N"], this->contentCurves["A/T/C/G/N"], this->contentCurves["GC"]
         * this->kmerMin, this->kmerMax, after summary, this->summarize will be set to true
         * @param forced forced to run summarize
         * this->summarize will only run if (this->summarized is false) || (this->summarized is true && force = false)
         */
        void summarize(bool forced = false);

        /** Generate json report
         * @param ofs std::ofstream to output report
         * @param padding padding in front of each line
         */
        void reportJson(std::ofstream& ofs, std::string padding);
        
        /** Generate Html report
         * @param ofs std::ofstream to output report
         * @param filteringType filtering type of report
         * @param readName library name 
         */
        void reportHtml(std::ofstream& ofs, std::string filteringType, std::string readName);
        

        /** Generate Quality part of the Html report
         * @param ofs std::ofstream to output report
         * @param filteringType filtering type of report
         * @param readName library name 
         */
        void reportHtmlQuality(std::ofstream& ofs, std::string filteringType, std::string readName);
        
        /** Generate Content part of the Html report
         * @param ofs std::ofstream to output report
         * @param filteringType filtering type of report
         * @param readName library name 
         */
        void reportHtmlContents(std::ofstream& ofs, std::string filteringType, std::string readName);
        
        /** Generate Kmer part of the Html report
         * @param ofs std::ofstream to output report
         * @param filteringType filtering type of report
         * @param readName library name 
         */
        void reportHtmlKmer(std::ofstream& ofs, std::string filteringType, std::string readName);
         
        /** Generate Over representation analysis  part of the Html report
         * @param ofs std::ofstream to output report
         * @param filteringType filtering type of report
         * @param readName library name 
         */
        void reportHtmlORA(std::ofstream& ofs, std::string filteringType, std::string readName);
        
        /** whether is long read fq
         * @return true if this->cycles > 300
         */
        bool isLongRead();
        
        /** get mean read length 
         * @return mean read length 
         */
        int getMeanLength();

        /** convert an array of value type T to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         */ 
        template<typename T>
        static std::string list2string(T* list, int size);
       
        /** convert an array of T values to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         * @param coords coordinates to define the ith value output
         * ith value output equals average of(list[coords[i-1]] ... list[coords[i]])
         */ 
        template<typename T>
        static std::string list2string(T* list, int size, size_t* coords);

    private:
        /** extend the array buffer for statistics longer
         * @param newBufLen the expected smallest buffer length
         * for performance, buffer length will increase to std::max(newBufLen + 100, 1.5 * newBufLen)
         */
        void extendBuffer(int newBufLen);
       
        /** make html popup value of kmer table
         * @param n the integer representation of kmer
         * @return popup value of kmer table
         */
        std::string makeKmerTD(size_t n);
        
        /** delete this->OverRepDist recources */ 
        void deleteOverRepSeqDist();

        /** test wheather count of seq is over represented
         * @param seq sequence to be evaluated
         * @param count seq count in sampling 
         * @return true if the seq is over represented
         */
        bool overRepPassed(const std::string& seq, size_t count);
};

#endif

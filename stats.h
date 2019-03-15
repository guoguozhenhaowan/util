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
        int kmerLen = 5;            ///< kmer length to calculate, default 5
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
        std::map<std::string, size_t> overRepSeq;     ///< map of <repSeq, repSeqCount>
        std::map<std::string, size_t*> overRepSeqDist;///< map of <repseq, dist*>

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
        void allocateRes();
        
        /** get cycle number of this->fqFile
         * @return cycle number of this->fqFile
         */
        int getCycles();
        
        /** get total raeds number of this->fqFile
         * @return total read number of this->fqFile
         */
        size_t getReads();
        
        /** get total base number of this->fqFile
         * @return total base number of this->fqFile
         */
        size_t getBases();
        
        /** get number of bases with quality equal or greater than 20 in this->fqFile
         * @return number of bases with quality equal or greater than 20 in this->fqFile
         */
        size_t getQ20();

        /** get number of bases with quality equal or greater than 30 in this->fqFile 
         * @return number of bases with quality equal or greater than 30 in this->fqFile
         */
        size_t getQ30();

        /** get number of base G and C in this->fqFile
         * @return number of base G and C in this->fqFile
         */
        size_t getGCNumber();

        /** do statistics of one read in this->fqFile
         * @param r pointer to object Read
         */
        void statRead(Read* r);

        /** merge a list of Stats objects into one
         * @param fqName fastqName of these Stats objects
         * @parma list a list of Stats objects
         * @return a merged Stats object
         */
        static Stats* merge(const std::string& fqName, std::vector<Stats*>& list);
        
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
        
        /** whether this->fqFile is long read fq
         * @return true if this->cycles > 300
         */
        bool isLongRead();
        
        /** get mean read length of this->fqFile
         * @return mean read length of this->fqFile 
         */
        int getMeanLength();

    public:
        /** convert an array of value type T to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         */ 
        template<typename T>
        static std::string list2string(T* list, int size);
       
        /** convert an array of T values to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         * @param coordinates to define the ith value output
         * ith value output equals average of(list[coords[i-1]] ... list[coords[i]])
         */ 
        template<typename T>
        static std::string list2string(T* list, int size, size_t* coords);

    public:
        
        /** extend the array buffer for statistics longer
         * @param newBufLen the expected smallest buffer length
         * for performance, buffer length will increase to std::max(newBufLen + 100, 1.5 * newBufLen)
         */
        void extendBuffer(int newBufLen);
       
        /** make html popup value of kmer table
         * @param i row index 
         * @param j col index
         * @return popup value of kmer table
         */
        std::string makeKmerTD(int i, int j);
        
        /** delete this->OverRepDist recources */ 
        void deleteOverRepSeqDist();

        /** test wheather count of seq is over represented
         * @param seq sequence to be evaluated
         * @count seq count in sampling of this->fqFile
         */
        bool overRepPassed(const std::string& seq, size_t count);
};

#endif

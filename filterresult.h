#ifndef FILTER_RESULT_H
#define FILTER_RESULT_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "common.h"
#include "stats.h"
#include "jsonutil.h"
#include "htmlutil.h"

/** filter options struct */
struct FilterOpt{
    bool filterShortRead;    ///< filter too short read if true
    bool filterLongRead;     ///< filter too long read if true
    bool trimAdapter;        ///< trim adapter sequnce if true
    bool baseCorrection;     ///< correct base if true
    bool complexityFilter;   ///< filter low complexity tread if true
    int minReadLen;          ///< read length lower than this will be treated as too short
    int maxReadLen;          ///< read length greater than this will be treated as too long
    FilterOpt() = default;
    ~FilterOpt() = default;
};

/** Class to store/calculate read filter results */
class FilterResult{
    public:
        FilterOpt *opt;                                         ///< pointer to FilterOpt object store various filter options
        bool paired;                                            ///< pairend filter results if true, else single end results
        size_t correctedReads;                                  ///< number of reads got bases corrected
        size_t correctedBases;                                  ///< number of bases got corrected
        size_t filterReadStats[compar::FILTER_RESULT_TYPES];    ///< reads filter status accumulation array
        size_t trimmedAdapterReads;                             ///< number of reads got adapter trimmed
        size_t trimmedAdapterBases;                             ///< number of bases got trimmed when do adapter trimming
        std::map<std::string, size_t> adapter1;                 ///< read1 adapter trimmed count map
        std::map<std::string, size_t> adapter2;                 ///< read2 adapter trimmed count map
        size_t* correctionMatrix;                               ///< 1x64 array to store various base to base correction accumulation numbers
        bool summarized;                                        ///< whether a FilterResult has been summarized(i.e correctedBases calculated)
    public:
        /** Construct FilterResult object to store/calculate filter results
         * @param opt FilterOpt object to store various filter arguments
         * @param paired if true the FilterResult object store pairend filter results
         */ 
        FilterResult(FilterOpt* opt, bool paired = false);
        
        /** Destruct FilterResult object, free correctionMatrix */
        ~FilterResult();

        /** get filterReadStats array(reads filter status accumulation array)
         * @return pointer to filterReadStats
         */
        inline size_t* getFilterReadStats(){
            return this->filterReadStats;
        }

        /** add filter status result to filter status result accumulation array
         * @param result filter status
         */
        void addFilterResult(int result);

        /** merge list of FilterResult object into one(combine filter results of various parts of a fastq file)
         * @param list a vector of pointer to FilterResult object
         * @return pointer to the merged FilterResult object
         */
        static FilterResult* merge(std::vector<FilterResult*>& list);
        
        /** Operator to output FilterResult to ostream
         * @param os reference of std::ostream object
         * @param re pointer of FilterResult object
         */
        friend std::ostream& operator<<(std::ostream& os, FilterResult* re);
        
        /** summary the FilterResult object, in case recalculate some parameters again and again
         * @param force if true, recalculate parameters
         */
        void summary(bool force = false);

        /** Update FilterResult when a adapter trimming happened
         * @param adapter adapter sequence trimmed
         * @param isR2 if true the adapter was trimmed from read2
         */
        void addAdapterTrimmed(const std::string& adapter, bool isR2 = false);

        /** Update FilterResult when adapter trimming happpened in pair end reads filtering
         * @param adapter1 adapter sequence trimmed from read1
         * @param adapter2 adapter sequence trimmed from read2
         */
        void addAdapterTrimmed(const std::string& adapter1, const std::string& adapter2);

        /** Report basic filter results to json file
         * @param ofs output file stream
         * @param padding padding before each json record
         */
        void reportJsonBasic(std::ofstream& ofs, const std::string& padding);

        /** Report basic filter results to html file
         * @param ofs output file stream
         * @param totalReads total reads of fastq
         * @param totalBases total bases of fastq
         */
        void reportHtmlBasic(std::ofstream& ofs, size_t totalReads, size_t totalBases);

        /** Report detailed trimmed adapter sequence and counts in json fomat
         * @param ofs output file stream
         * @param adapterCounts adapter sequence counts map
         */
        void reportAdaptersJsonDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts);
        
        /** Report summary trimmed adapter sequence and counts in json format
         * @param ofs output file stream
         * @param padding padding before each json record
         * @param adapterSeq1 read1 adapter sequence
         * @param adapterSeq2 read2 adapter seqeunce
         */
        void reportAdaptersJsonSummary(std::ofstream& ofs, const std::string& padding, const std::string& adapterSeq1, const std::string& adapterSeq2="");
        
        /** Report detailed trimmed adapter sequence and counts in html format
         * @param ofs output file stream
         * @param adapterCounts adapter sequence counts map
         * @param totalBases totalBases total bases of fastq
         */
        void reportAdaptersHtmlDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts, size_t totalBases);

        /** Report summary trimmed adapter sequence and counts in html format
         * @param ofs output file stream
         * @param totalBases totalBases total bases of fastq
         */
        void reportAdaptersHtmlSummary(std::ofstream& ofs, size_t totalBases);

        /** get correction matrix(1x64 array)
         * @return pointer to correctionMatrix
         */
        size_t* getCorrectionMatrix(){
            return this->correctionMatrix;
        }

        /** get number of bases corrected
         * @return number of based corrected
         */
        size_t getTotalCorrectedBases();

        /** add one case of base correction to FilterResult object, \\nupdate the count in correctionMatrix,\\n
         * correctionMatrix[(from & 0x07) * 8 + (to & 0x07)] += 1
         * @param from source uncorrected base
         * @param to destination corrected base
         */
        void addCorrection(char from, char to);
       
        /** get base correction numbers in the pattern (from->to)
         * @param from source uncorrected base
         * @param to destination corrected base
         * @return base correction numbers in the pattern(from->to), just correctionMatrix[(from & 0x07) * 8 + (to & 0x07)]
         */ 
        size_t getCorrectionNum(char from, char to);

        /** increase number of corrected reads
         * @param count number of corrected reads to increase
         */ 
        void incCorrectedReads(int count);
};

#endif

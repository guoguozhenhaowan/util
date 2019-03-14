#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <map>
#include <set>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "read.h"
#include "util.h"
#include "fqreader.h"

/** class to hold various functions to evaluate sequence information */
class Evaluator{
    public:
        std::string r1File;  ///< read1 filename
        std::string r2File;  ///< read2 filename
    public:
        /** Construct a evaluator with read1/2 filename
         * @param r1f read1 filename
         * @param r2f read2 filename
         */
        Evaluator(const std::string& r1f, const std::string& r2f = "") : r1File(r1f), r2File(r2f){}
        ~Evaluator(){}

        void evaluateReadNum(size_t& readNum);
        std::string evalAdapterAndReadNumDepreciated(size_t& readNum);
        std::string evalAdapterAndReadNum(size_t& readNum, bool isR2);
        
        /** Test whether the read is sequenced from a TwoColorSystem machine
         * @return true if read name starts with "\@NS", "\@NB" or "\@A0"
         * which corresponding to NEXTSEQ500, NEXTSEQ550 and NOVASEQ machine
         */
        bool isTwoColorSystem();
        
        /** Evaluate the read1/2 length if exists
         * @param r1Len value to store read1 length
         * @param r2Len value to store read2 length
         */
        void evaluateSeqLen(int& r1Len, int& r2Len);
        
        /** Evaluate the fastq file over represented sequences based on at most 151 * 10000 bases
         * count subsequences of length in {10, 20, 40, 100, min(150, seqLen -2)} in each reads(total bases less than 151 * 10000)
         * a subsequence will be considered as over represented if 
         * (length >= seqLen - 1 && count >= 3) || (length >= 100 && count >= 5) || 
         * (length >= 40 && count >= 20) || (length >= 20 && count >= 100 || (length >= 10 && count >= 500)
         * remove substrings in the map if the count of substring is less than 10 * count of the string contains it
         * @param filename fastq filename
         * @param hotSeqs map to store over represented sequences
         * @param seqLen max sequence length of fastq filename
         */ 
        void computeOverRepSeq(const std::string& filename, std::map<std::string, size_t>& hotSeqs, const int& seqLen);
        
        /** Evaluate the max read length of a fastq file
         * @param filename fastq file
         * @return the max read length of this fastq file based on the first 1000 reads length
         */
        int computeSeqLen(const std::string& filename);

        static std::string matchKnownAdapter(std::string& seq);
    
    public:
        std::string int2seq(unsigned int val, int seqLen);
        int seq2int(std::string& seq, int pos, int seqLen, int lastVal = -1);
        std::string getAdapterWithSeed(int seed, Read** loadedReads, size_t records, int keyLen);
};

#endif

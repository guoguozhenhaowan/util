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
#include "knownadapters.h"
#include "nucleotidetree.h"

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
        
        /** Destroy a evaluator */
        ~Evaluator(){}

        /** evaluate number of reads in fastq this->r1File 
         * based on at most 512 * 1024 reads and at most 151 * 512 * 1024 bytes read
         * @param readNum store the estimated reads number 
         */
        void evaluateReadNum(size_t& readNum);
        
        /** evaluate the total reads number of this->r1File or this->r2File 
         * and the adapter sequences based on at most 256*1024 reads or 151 * 256 * 1024 bytes read in
         * to estimate adapter sequence, there must be at least 10000 valid records read in
         * @param readNum store the estimated reads number
         * @param isR2 evaluate this->r2File if true
         * @param trim cycles to trim from end(3') before evaluate adaptors
         * @return estimated adapter sequences if successful
         */
        std::string evalAdapterAndReadNum(size_t& readNum, bool isR2, size_t trim);
        
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

        /** Match a sequence seq agains known illumina adapters stored in knownadapters.h
         * seq length must be equal or greater than the matched adapter 
         * and match happens from index 0 between seq and adapter exactly to the end of adapter
         * not any mismatch allowes 
         * @param seq a nucleotide sequence(may generated from NeucleotideTree)
         * @return matched adapter sequence if successfuly or "" if failed
         */
        static std::string matchKnownAdapter(std::string& seq);
    
    public:
        /** Convert a size_t value back into a string consists of ATCG with length seqLen
         * @param val the encoded value of a nucleotide sequence
         * @param seqLen the encoded nucleotide sequence length
         * @return the nucleotide sequence of val encoded
         */
        std::string int2seq(size_t val, int seqLen);
        
        /** convert the substring started at pos to pos + seqLen - 1 of a string to int
         * each 2 bits represents a nucleotide[A->00, T->01, C->10, G->11]
         * from lowest bits to highest bits represent seq[pos + seqLen - 1] to seq[pos]
         * if nucleotide other than ATCG occurs in this substring, -1 will return 
         * this sliding window algorithm will have a complexity of o(n)
         * @param seq the whole string
         * @param pos ini position of substring of seq
         * @param seqLen length of substring to be converted
         * @param lastVal the int representation of  substring started at pos - 1 to pos + seqLen - 2 of string
         * @return int representation of the substring started at pos to pos + seqLen - 1 of string seq, may be -1
         */
        int seq2int(std::string& seq, int pos, int seqLen, int lastVal = -1);

        /** get adapter sequence from a seed sequence represented by keyLen bits int
         * a seed sequence is a substring of read(trimed with trim) with length keyLen
         * and its in the top10 counts(excluding low complexity, high GC, GGGG** seq, count < 10 
         * or count < n * p ( n = FOLD_THRESHOLD, p = total/size) 
         * @param seed int representation of the seed sequence
         * @param loadedReads loaded Reads to get the seed 
         * @param records number of Reads loaded
         * @param keyLen length of subsquence
         * @param trim trim length of read from 3' end
         * @return adapter sequence detected
         */
        std::string getAdapterWithSeed(int seed, Read** loadedReads, size_t records, int keyLen, int trim);
};

#endif

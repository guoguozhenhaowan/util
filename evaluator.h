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
        std::string fqFile;   ///< fastq readfilename
        int trimLen;          ///< trim length from 3' end before evaluation
        int readLen;          ///< maximum read length
        size_t readNum;       ///< estimated read numbers
        bool twoColorSystem;  ///< estimated this fastq come frm a two color sequencer
        std::string adapter;  ///< estimated adapter sequence
        bool illuminaAdapter; ///< is the adapter estimated is a standard known illumina adapter
    public:
        /** Construct a evaluator with read1/2 filename
         * @param trlen trim length from 3' end before evaluation
         * @param filename fastq filename
         */
        Evaluator(const std::string& filename, const int& trlen) : trimLen(trlen), fqFile(filename){
            this->readLen = 0;
            this->readNum = 0;
            this->twoColorSystem = false;
            this->illuminaAdapter = false;
            this->adapter = "";
            this->evaluateReadLen();
            this->evaluateReadNum();
            this->evaluateTwoColorSystem();
        }
        
        /** Destroy a evaluator */
        ~Evaluator(){}

        /** get the estimated read number of this fastq file
         * @return the estimated read number of this fastq file
         */
        inline size_t getReadNum(){
            return this->readNum;
        }
        
        /** get the estimated adapter sequence of this fastq file
         * @return the estimated adapter sequence of this fastq file
         */
        inline std::string getAdapter(){
            this->evaluateAdapterSeq();
            return this->adapter;
        }

        /** whether the adapter estimated is illumina know adapter
         * @return true if the adapter estimated is illumina know adapter
         */
        inline bool isIlluminaAdapter(){
            return this->illuminaAdapter;
        }

        /** Test whether the read is sequenced from a TwoColorSystem machine
         * @return true if the read is sequenced from a TwoColorSystem machine
         */
        inline bool isTwoColorSystem(){
            return this->twoColorSystem;
        }
        
        /** Get the estimated read length of the fastq file
         * @return the estimated read length of the fastq file
         */
        inline int getReadLen(){
           return this->readLen;
        }

        /** Evaluate the fastq file over represented sequences based on at most this->readLen * 10000 bases
         * count subsequences of length in {10, 20, 40, 100, min(150, this->readLen -2)} in each reads(total bases less than this->readLen * 10000)
         * a subsequence will be considered as over represented if 
         * (length >= this->readLen - 1 && count >= 3) || (length >= 100 && count >= 5) || 
         * (length >= 40 && count >= 20) || (length >= 20 && count >= 100 || (length >= 10 && count >= 500)
         * remove substrings in the map if the count of substring is less than 10 * count of the string contains it
         * @param hotSeqs map to store over represented sequences
         */ 
        void computeOverRepSeq(std::map<std::string, size_t>& hotSeqs);

        /** Match a sequence seq agains known illumina adapters stored in knownadapters.h
         * seq length must be equal or greater than the matched adapter 
         * and match happens from index 0 between seq and adapter exactly to the end of adapter
         * not any mismatch allowes 
         * @param seq a nucleotide sequence(may generated from NeucleotideTree)
         * @return matched adapter sequence if successfuly or "" if failed
         */
        static std::string matchKnownAdapter(const std::string& seq);

        /** Convert a size_t value back into a string consists of ATCG with length seqLen
         * @param val the encoded value of a nucleotide sequence
         * @param seqLen the encoded nucleotide sequence length
         * @return the nucleotide sequence of val encoded
         */
        static std::string int2seq(size_t val, int seqLen);
        
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
        static int seq2int(const std::string& seq, int pos, int seqLen, int lastVal = -1);

    private:
        /** evaluate number of reads in fastq file
         * based on at most 512 * 1024 reads and at most this->readLen * 512 * 1024 bytes read
         */
        void evaluateReadNum();

        /** evaluate the adapter sequences based on at most 256*1024 reads or this->readLen * 256 * 1024 bytes read in
         * to estimate adapter sequence, there must be at least 10000 valid records read in
         */
        void evaluateAdapterSeq();

        /** Test whether the read is sequenced from a TwoColorSystem machine
         * which corresponding to NEXTSEQ500, NEXTSEQ550 and NOVASEQ machine
         */
        void evaluateTwoColorSystem();
        
        /** Evaluate the max read length of a fastq file */
        void evaluateReadLen();
        
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
        std::string getAdapterWithSeed(int seed, Read** loadedReads, long records, int keyLen, int trim);
};

#endif

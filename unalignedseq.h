#ifndef UNALIGNED_SEQ_H
#define UNALIGNED_SEQ_H

#include <cstring>
#include <vector>
#include <iostream>

/** Structure to hold unaligned sequence (name and bases) */
class UnalignedSeq{
    public:
        std::string mName;    ///< name of the contig
        std::string mComment; ///< comment of the contig
        std::string mSeq;     ///< sequence of the contig (upper-case ACTGN)
        std::string mQual;    ///< quality scores
        char mStrand;         ///< strand of the sequence. Default is '*'
    
    public:
        /** Construct an empty sequence */
        UnalignedSeq(){}
       
        /** Construct an unaligned sequence with name and sequence
         * @param n name of the sequence 
         * @param s sequence, stored as ACTG or N characters
         */
        UnalignedSeq(const std::string& n, const std::string& s) : mName(n), mComment(std::string()), mSeq(s), mQual(std::string()), mStrand('*') {}

        /** Construct an unaligned mSequence with name, sequence and quality score
         * @param n name of the sequence 
         * @param s sequence, stored as ACTG or N characters
         * @param q quality string
         */
        UnalignedSeq(const std::string& n, const std::string& s, const std::string& q) : mName(n), mComment(std::string()),  mSeq(s), mQual(q), mStrand('*') {}

        /** Construct an unaligned sequence with name, sequence, quality score and strand
         * @param n name of the sequence 
         * @param s sequence, stored as ACTG or N characters
         * @param q quality string
         * @param t strand of the sequence, one of '*', '+', '-'
         */
        UnalignedSeq(const std::string& n, const std::string& s, const std::string& q, char t) : mName(n), mComment(std::string()), mSeq(s), mQual(q), mStrand(t) {}
        
        /** Output an unaligned sequence to ostream
        * @param os ostream
        * @param us UnalignedSeq
        */
        friend std::ostream& operator<<(std::ostream& os, const UnalignedSeq& us){
            os << "@" << us.mName << " " << us.mComment << "\n";
            os << us.mSeq << "\n";
            os << us.mStrand << "\n";;
            os << us.mQual << "\n";
            return os;
        }
};

#endif

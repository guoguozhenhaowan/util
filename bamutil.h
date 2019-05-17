#ifndef BAMUTIL_H
#define BAMUTIL_H

#include <sstream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <string>
#include <cstring>
#include <algorithm>
#include "util.h"
#include "htslib/sam.h"

/** some usefule functions to operate bam file */
namespace bamutil{
    /** get read name of an alignment record
     * @param b pointer to bam1_t struct
     * @return read name of the alignment
     */ 
    std::string getQName(bam1_t* b);
    
    /** get barcode sequence of an alignment record\n
     * the barcode sequence should be put into fastq read header as BC:Z:BarcodeSeq\n
     * and then use aligner to put this OX:Z tag into bam record such as bwa mem -C\n
     * @param b pointer to bam1_t struct
     * @return barcode sequence of the alignment if OX tag parsed properly else empty string
     */
    std::string getBarcode(bam1_t* b);
    
    /** get read sequence of an alignment record
     * @param b pointer to bam1_t struct
     * @return read seqence of the alignment
     */
    std::string getSeq(bam1_t* b);
    
    /** get quality string of an alignment record
     * @param b pointer to bam1_t struct
     * @return quality string of the alignment
     */
    std::string getQual(bam1_t* b);
    
    /** get cigar string of an alignment record
     * @param b pointer to bam1_t struct
     * @return cigar string
     */
    std::string getCigar(bam1_t* b);
    
    /** output an alignment record to std::cerr
     * @param b pointer to bam1_t struct
     */
    void dump(bam1_t* b);

    /** test wheather an alignment record is part of another alignment record
     * @param part pointer to bam1_t struct
     * @param whole pointer to bam1_t struct
     * @param isLeft compare from left to right if true,\n
     * that is to say part and whole are mapped to reference at the same start position
     */
    bool isPartOf(bam1_t* part, bam1_t* whole, bool isLeft);
    
    /** output bam header into std::cerr
     * @param hdr pointer to struct bam_hdr_t
     */
    void dumpHeader(bam_hdr_t* hdr);
    
    /** return length of reference consumed in range [b->core.pos to b->core.pos + offset)
     * @param b pointer to bam1_t struct
     * @param offset offset from read mapping starting position of b
     * @return length of reference consumed until [b->core.pos to b->core.pos + offset)\n
     * or -1 if offset is invalid or cigar is invalid(unmapped)
     */
    int getRefOffset(bam1_t* b, int offset);
    
    /** copy read name of one alignment record to another
     * @param from pointer to bam1_t struct
     * @param to pointer to bam1_t struct
     */
    void copyQName(bam1_t* from, bam1_t* to);

    /** modify read name of one alignment record 
     * @param b pointer to bam1_t struct
     * @param n new name to be used in this read
     */
    void changeQname(bam1_t* b, const std::string n);

    /** modify the ith base of one alignment record
     * @param b pointer to bam1_t struct
     * @param i position at which to change nucleotide(0 based)
     * @param ch nucleotide to change to
     */
    void changeSeq(bam1_t* b, int i, char ch);
    

    /** modify the ith base of one alignment record
     * @param b pointer to bam1_t struct
     * @param i position at which to change nucleotide(0 based)
     * @param ch nucleotide to change to
     */
    void changeSeq(bam1_t* b, int i, uint8_t ch);
    
    /** test thether the alignment is primary alignment
     * @return true if this alignment record is primary
     */
    bool isPrimary(bam1_t* b);
    
    /** test whether read1/2 in this alignment mapped properly
     * @parm b pointer to bam1_t struct
     */
    bool isProperPair(bam1_t* b);
    
    /** calculate the rightmost base position of an alignment on the reference genome 
     * @param b pointer to bam1_t struct
     * @return rthe coordinate of the first base after the alignment, 0-based\n
     * or 0 if unmapped
    */
    int getRightRefPos(bam1_t* b);
    
    /** get the first matched status
     * @param b pointer to bam1_t struct
     * @MatchOffset length of sequence consumed in read befor the first cigar M met
     * @param MatchLen the first cigar M consumed length
     */
    void getFirstMatchOffsetAndLen(bam1_t* b, int& MatchOffset, int& MatchLen);
    
    /** convert a 4 bits in bam alignment record data sequence array into a base\n
     * this sequence can be got by calling bam_get_seq(b)\n
     * @param val a 4 bits value
     * @return corresponding nucleotide 1->A, 2->C, 4->G, 8->T, 15->N, otherwise(N)
     */
    char fourbits2base(uint8_t val);

    /** convert a base character to a four bits uint8_t
     * @param chr a nucleotide character
     * @return a uint8_t value which can be used in bam alignment record data sequence array to represent this nucleotide\n
     * A->1, C->2, G->4, T->8, N->15
     */
    uint8_t base2fourbits(char chr);
    
    /** get edit distance of this alignment record
     * @param b pointer to bam1_t struct
     * @return edit distance of this alignment record if NM tag parsed properly, else 0
     */
    int getED(bam1_t* b);

    /** get quality sum of a read
     * @param b pointer to bam1_t struct
     * @return quality sum of read in b
     */
    uint8_t getReadQualSum(bam1_t* b);

    /** get effective bases
     * @param b pointer to bam1_t struct
     * @return mapped bases number
     */
    uint8_t getMappedBases(bam1_t* b);

    /** get mismatches in range
     * @param b pointer to bam1_t struct
     * @return mismatches positions on reference(0 based)
     */
    std::vector<int32_t> getMismatchPos(bam1_t* b);

    /** a bam record is paired
     * @param b pointer to bam1_t struct
     * @return true if b is paired
     */
    inline bool bamIsPaired(const bam1_t* b){
        return b->core.flag & BAM_FPAIRED;
    }

    /** a bam record is proper paired 
     * @param b pointer to bam1_t struct
     * @return true if b is properly paired
     */
    inline bool bamIsProper(const bam1_t* b){
        return b->core.flag & BAM_FPROPER_PAIR;
    }

    /** a bam record is mapped
     * @param b pointer to bam1_t struct
     * @return true if b is mapped
     */
    inline bool bamIsMapped(const bam1_t* b){
        return !(b->core.flag & BAM_FUNMAP);
    }

    /** a bam record is unmapped
     * @param b pointer to bam1_t struct
     * @return true if b is unmapped
     */
    inline bool bamIsUnmapped(const bam1_t* b){
        return b->core.flag & BAM_FUNMAP;
    }

    /** a bam record is for read1
     * @param b pointer to bam1_t struct
     * @return true if b is alignment record of read1
     */
    inline bool bamIsRead1(const bam1_t* b){
        return b->core.flag & BAM_FREAD1;
    }

    /** a bam record is for read2
     * @param b pointer to bam1_t struct
     * @return true if b is alignment record of read2
     */
    inline bool bamIsRead2(const bam1_t* b){
        return b->core.flag & BAM_FREAD2;
    }

    /** get indel position on read of an alignment record
     * @param p pointer to bam1_t struct
     * @param indelPosVec list of indel positions on read, 0 based
     */
    inline void getIndelPos(const bam1_t* b, std::vector<int>& indelPosVec);

    /** check if indels are around the current position
     * @param pos position on read, 0-based
     * @param range range to be checked around pos
     * @param ivec indel position on read of a alignment record
     * @return true if a indel occurs around pos±range inclusive[]
     */
    inline bool posIsAdjacentIndel(int pos, int range, const std::vector<int>& ivec);

    /** get compacted mismatches of a read against the reference in an alignment record
     * @param b pointer to bam1_t struct
     * @return compacted mismatches of b(each indel is counted as one mismatch, ambiguous N mismatch discounted)
     */
    inline int getCompactMismatch(bam1_t* b);

    /** get max read length in a bam by iterate at most the first 10000 record
     * @param bamFile bam filename
     * @return max read length scanned
     */
    int getBamReadLength(const std::string& bamFile);

    /** check if a alignment record in a bam_pileup1_t pileup near some range of read 3' end
     * @param p pointer to bam_pileup1_t struct
     * @param range range to be checked near 3'end of read
     * @return true if pipeup position is in read 3'end-range, inclusive[]
     */
    bool readRearPartPileup(const bam_pileup1_t* p, int range);

    /** get average base quality around predefined range of of pileup posision
     * @pram p pointer to bam_pileup1_t 
     * @param range max number of bases to counted around p->qpos
     * @return average base quality around qpos±range
     */
    float getAverageBaseQualityAroundPileupPos(const bam_pileup1_t* p, int range);

    /** get at most insertion sequence at pileup position of a bam record
     * @param p pointer to bam_pileup1_t 
     * @param maxLen max insertion sequence to extract
     * @return insertion sequence extracted
     */
    std::string getInsertionSeq(const bam_pileup1_t* p, int maxLen);

    /** get softclip length of left/right softclip
     * @param b pointer to bam1_t
     * @return <leftClipLength, rightClipLength>
     */
    std::pair<int, int> getSoftClipLength(bam1_t* b);

    /** get minimal distance from pileup position to the end of read
     * @param p pointer to bam_pileup1_t
     * @return minimal distance from p->qpos to ends of read
     */
    int getMinDistToReadEnd(const bam_pileup1_t* p);
}

#endif

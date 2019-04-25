#ifndef ONLINEBWA_H
#define ONLINEBWA_H

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include <cassert>
#include "util.h"
#include "bwa/bwa.h"
#include "bwa/bwt.h"
#include "bwa/kseq.h"
#include "bwa/utils.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "unalignedseq.h"

extern "C"{
    int is_bwt(ubyte_t *T, int n);
    KSEQ_DECLARE(gzFile)
}

class OnlineBWA{
    public:
        mem_opt_t* mMemOpt; ///< pointer to memory options for bwa alignment
        bool mCopyComment;  ///< copy fasta/q comment to SAM output if true
        bwaidx_t* mIndex;   ///< pointer to bwa index structure
        
        /** OnlineBWA constructor */ 
        OnlineBWA();
        
        /** OnlineBWA destructor */
        ~OnlineBWA();

    public:
        /** convert a bns to string
         * @return pointer to bam_hdr_t
         */
        bam_hdr_t* getBamHeader();
        
        /** convert pac tos bwt_t
         * @param pac uint8_t array
         * @param bwt_seq_lenr bwt sequence length
         * @return pointer to bwt_t
         */
        bwt_t* pac2bwt(const uint8_t* pac, int bwt_seq_lenr);
        
        /** add an annotation sequence
         * @param name reference of annotation name
         * @param seq reference of annotation sequence
         * @param ann pointer to bntann1_t struct
         * @param offset offset of bntann1_t to add annotation
         */
        void addAnns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset);

        /** bwa mem add1 function */
        uint8_t* addSeq2pac(const kseq_t* seq, bntseq_t* bns, uint8_t* pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q);

        /** make the pac structure for a bunch of reads
         * @param v vector of UnalignedSeq 
         * @param addReverse add reverse sequence to pac if true
         */
        uint8_t* makePac(std::vector<UnalignedSeq> v, bool addReverse);

        /** write pac part of the index to file
         * @param file filename to write pac to
         */
        void writePacToFile(const std::string& file) const;

        /** align a sequence to reference
         * @param seq nucleotide sequence to do align
         * @name sequence name
         * @result alignment result vector
         */
        void alignSeq(const std::string& seq, const std::string& name, std::vector<bam1_t*>& result);

        /** align a UnalignedSeq to reference
         * @param us reference of UnalignedSeq object
         * @param result alignment result vector
         */
        void alignSeq(const UnalignedSeq& us, std::vector<bam1_t*>& result);

        /** align a kseq_t read to reference
         * @param seq pointer to kseq_t struct
         * @param result alignment result vector
         */
        void alignSeq(const kseq_t* seq, std::vector<bam1_t*>& result);

        /** construct index from a list of UnalignedSeq
         * @param uv vector of UnalignedSeq
         */
        void constructIndex(const std::vector<UnalignedSeq>& uv);

        /** load external bwt index
         * @param file index file
         */
        void loadIndex(const std::string& file);

        /** write index to file
         * @param file output file of index
         */
        void writeIndex(const std::string& file);

        /** set the gap open penalty
         * @param gapOpenPenalty gap open penalty, default 6
         */
        inline void setGapOpenPenalty(int gapOpenPenalty){
            assert(gapOpenPenalty > 0);
            mMemOpt->o_del = mMemOpt->o_ins = gapOpenPenalty;
        }

        /** set the gap extension penalty
         * @param gapExtPenalty gap extension penalty, default 1
         */
        inline void setGapExtendPenalty(int gapExtPenalty){
            assert(gapExtPenalty > 0);
            mMemOpt->e_del = mMemOpt->e_ins = gapExtPenalty;
        }

        /** set mismatch penalty
         * @param mismatchPenaly mismatch penalty, default 4
         */
        inline void setMismatchPenalty(int mismatchPenaly){
            assert(mismatchPenaly > 0);
            mMemOpt->b = mismatchPenaly;
            bwa_fill_scmat(mMemOpt->a, mMemOpt->b, mMemOpt->mat);
        }

        /** set the reseed trigger
         * @param reseed look for internal seeds inside a seed longer than seedlength * reseed, default 1.5
         */
        inline void setReseedTriger(float reseed){
            assert(reseed > 0);
            mMemOpt->split_factor = reseed;
        }

        /** set SW alignment bandwidth
         * @param width SW alignment bandwidth, default 100
         */
        inline void setBandWidth(int width){
            assert(width > 0);
            mMemOpt->w = width;
        }

        /** set the SW alignment off-diagonal X-dropoff
         * @param dropOff off-diagonal X-dropoff, default 100
         */
        inline void setXDropoff(int dropOff){
            assert(dropOff > 0);
            mMemOpt->zdrop = dropOff;
        }

        /** set the 3' clipping penalty
         * @param p3ClipPenalty penalty for 3'-end clipping, default 5
         */
        inline void set3PrimeClipPenalty(int p3ClipPenalty){
            assert(p3ClipPenalty > 0);
            mMemOpt->pen_clip3 = p3ClipPenalty;
        }

        /** set the 5' clipping penalty
         * @param p5ClipPenalty penalty for 5'-end clipping, default 5
         */
        inline void set5PrimeClipPenalty(int p5ClipPenalty){
            assert(p5ClipPenalty > 0);
            mMemOpt->pen_clip5 = p5ClipPenalty;
        }

       /** set the match score, this should be set first as it will scale penalty options 
        * @param matchScore score for a sequence match, which scales options -TdBOELU unless overridden
        */
        inline void setMatchScore(int matchScore){
            assert(matchScore > 0);
            mMemOpt->b *= matchScore;
            mMemOpt->T *= matchScore;
            mMemOpt->o_del *= matchScore;
            mMemOpt->o_ins *= matchScore;
            mMemOpt->e_del *= matchScore;
            mMemOpt->e_ins *= matchScore;
            mMemOpt->zdrop *= matchScore;
            mMemOpt->pen_clip3 *= matchScore;
            mMemOpt->pen_clip5 *= matchScore;
            mMemOpt->pen_unpaired *= matchScore;
            mMemOpt->a = matchScore;
        }
};

#endif

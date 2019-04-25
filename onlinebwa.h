#ifndef ONLINEBWA_H
#define ONLINEBWA_H

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include "util.h"
#include "bwa/bwa.h"
#include "bwa/bwt.h"
#include "bwa/kseq.h"
#include "bwa/utils.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "htslib/sam.h"
#include "unalignedseq.h"

int is_bwt(ubyte_t *T, int n);
KSEQ_DECLARE(gzFile)

class OnlineBWA{
    public:
        mem_opt_t* mMemOpt; ///< pointer to memory options for bwa alignment
        bool mCopyComment; ///< copy fasta/q comment to SAM output if true
        bwaidx_t* mIndex; ///< pointer to bwa index structure
        
        OnlineBWA();
        ~OnlineBWA();

    public:
        /** convert a bns to string
         * @param bns pointer to bntseq_t struct
         * @return string representation of bns
         */
        std::string getSamHeader(const bntseq_t* bns);
        
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
};

#endif

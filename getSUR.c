#include <htslib/sam.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char** argv){
    if(argc < 3){
        printf("%s <in.bam> <out.bam>\n", argv[0]);
        return 0;
    }
    samFile* ofp = sam_open(argv[2], "w");
    samFile* ifp = sam_open(argv[1], "r");
    bam_hdr_t* hdr = sam_hdr_read(ifp);
    assert(sam_hdr_write(ofp, hdr) >= 0);
    bam1_t* b = bam_init1();
    uint16_t SUR_MASK = (BAM_FUNMAP & (~BAM_FMUNMAP));
    while(sam_read1(ifp, hdr, b) >= 0){
        if(b->core.flag & SUR_MASK){
            assert(sam_write1(ofp, hdr, b) >= 0);
        }
    }
    sam_close(ofp);
    sam_close(ifp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
}

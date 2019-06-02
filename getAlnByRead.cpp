#include <htslib/sam.h>
#include <iostream>
#include <string>

int main(int argc, char** argv){
    if(argc <= 3){
        std::cout << argv[0] << " <inbam> <readname> <outbam>" << std::endl;
        return 0;
    }
    samFile* ifp = sam_open(argv[1], "r");
    bam_hdr_t* h = sam_hdr_read(ifp);
    bam1_t* b = bam_init1();
    samFile* ofp = sam_open(argv[3], "wb");
    assert(sam_hdr_write(ofp, h) >= 0);
    while(sam_read1(ifp, h, b) >= 0){
        if(std::strcmp(bam_get_qname(b), argv[2]) == 0){
            assert(sam_write1(ofp, h, b) >= 0);
        }
    }
    sam_close(ofp);
    sam_close(ifp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

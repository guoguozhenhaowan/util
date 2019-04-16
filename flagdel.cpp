#include "htslib/sam.h"
#include <iostream>
#include <cassert>

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << argv[0] << " <inbam> <outbam> <flag>" << std::endl;
        return 0;
    }
    char* inbam = argv[1];
    char* outbam = argv[2];
    uint16_t flag = std::atoi(argv[3]);
    std::cout << "flag to be deleted: " << flag << std::endl;
    uint16_t mask = ~flag;

    samFile* ifh = sam_open(inbam, "r");
    bam_hdr_t* h = sam_hdr_read(ifh);
    bam1_t* b = bam_init1();
    
    samFile* ofh = sam_open(outbam, "wb");
    assert(sam_hdr_write(ofh, h) >= 0);

    uint32_t records = 0;
    while(sam_read1(ifh, h, b) >= 0){
        if(b->core.flag & flag){
            ++records;
            b->core.flag &= mask;
        }
        assert(sam_write1(ofh, h, b) >= 0);
    }
    sam_close(ifh);
    sam_close(ofh);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    std::cout << "flag deleted records: " << records << std::endl;
}


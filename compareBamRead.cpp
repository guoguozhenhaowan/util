#include "htslib/sam.h"
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << argv[0] << " <bam1> <bam2> " << std::endl;
        return 0;
    }
    char* inBam1 = argv[1];
    char* inBam2 = argv[2];
    samFile* fp = sam_open(inBam1, "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    std::map<std::string, std::pair<int, int>> bamReads;
    bam1_t* b = bam_init1();
    uint16_t mask = (BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    while(sam_read1(fp, h, b) >= 0){
        if(b->core.flag & mask){
            continue;
        }
        char* qname = bam_get_qname(b);
        if(bamReads.find(qname) == bamReads.end()){
            bamReads[qname] = {0, 0};
        }
        if(b->core.flag & BAM_FREAD1){
            bamReads[qname].first = 1;
        }
        if(b->core.flag & BAM_FREAD2){
            bamReads[qname].second = 1;
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    std::stringstream sReadOtherDiff;;
    std::stringstream sReadInBam2Only;
    fp = sam_open(inBam2, "r");
    h = sam_hdr_read(fp);
    while(sam_read1(fp, h, b) >= 0){
        if(b->core.flag & mask){
            continue;
        }
        char* qname = bam_get_qname(b);
        if(bamReads.find(qname) == bamReads.end()){
            sReadInBam2Only << qname << "\n";
            continue;
        }
        if(b->core.flag & BAM_FREAD1){
            bamReads[qname].first += 1;
        }
        if(b->core.flag & BAM_FREAD2){
            bamReads[qname].second += 1;
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    for(auto& e: bamReads){
        if(e.second.first + e.second.second < 4){
            sReadOtherDiff << e.first << "\n";
        }
    }
    std::cout << "Reads Only in " << inBam2 << ":\n";
    std::cout << sReadInBam2Only.str() << std::endl;
    std::cout << "Reads Diff Other: \n";
    std::cout << sReadOtherDiff.str() << std::endl;
}

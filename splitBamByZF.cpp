#include "htslib/sam.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << argv[0] << " <inbam> <outdir> " << std::endl;
        return 0;
    }
    
    char* inBam = argv[1];
    std::string outDir = argv[2];
    std::map<std::string, std::vector<bam1_t*>> recMap; 

    const char* tag = "ZF";
    samFile* ifp = sam_open(inBam, "r");
    bam_hdr_t* ibh = sam_hdr_read(ifp);
    bam1_t* rec = bam_init1();
    while(sam_read1(ifp, ibh, rec) >= 0){
        uint8_t* zftab = bam_aux_get(rec, tag);
        char* zfstr = bam_aux2Z(zftab);
        if(recMap.find(zfstr) == recMap.end()){
            recMap[zfstr] = std::vector<bam1_t*>();
            recMap[zfstr].push_back(rec);
        }else{
            recMap[zfstr].push_back(rec);
        }
        rec = bam_init1();
    }
    sam_close(ifp);
    for(auto& e: recMap){
        std::string outFile = outDir + "/" + e.first + ".bam";
        samFile* ofp = sam_open(outFile.c_str(), "wb");
        assert(sam_hdr_write(ofp, ibh) >= 0);
        for(auto& f: e.second){
            assert(sam_write1(ofp, ibh, f) >= 0);
            bam_destroy1(f);
            f = NULL;
        }
        sam_close(ofp);
        assert(sam_index_build(outFile.c_str(), 0) >= 0);
    }
    bam_hdr_destroy(ibh);
}

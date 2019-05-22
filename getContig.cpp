#include "htslib/faidx.h"
#include <iostream>

int main(int argc, char** argv){
    if(argc == 1){
        std::cout << argv[0] << " <ref> <contig> <start> <len>" << std::endl;
        return 0;
    }
    char* refFile = argv[1];
    char* contig = argv[2];
    int start = std::atoi(argv[3]);
    int len = std::atoi(argv[4]);
    faidx_t* fai = fai_load(refFile);
    int glen = 0;
    char* seq = faidx_fetch_seq(fai, contig, start, start + len - 1, &glen);
    std::cout << ">" << contig << "\n" << seq << std::endl;
    fai_destroy(fai);
}

#include "htslib/faidx.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << argv[0] << " <ref> <contig> [[start] [len]]" << std::endl;
        return 0;
    }
    char* refFile = argv[1];
    faidx_t* fai = fai_load(refFile);
    assert(fai);
    char* contig = argv[2];
    int start = 0;
    if(argc > 3) start = std::atoi(argv[3]);
    int len = faidx_seq_len(fai, contig);
    if(argc > 4) len = std::atoi(argv[4]);
    int glen = 0;
    char* seq = faidx_fetch_seq(fai, contig, start, start + len - 1, &glen);
    std::cout << "seq name quired: " << contig << std::endl;
    std::cout << "seq length got : " << glen << std::endl;
    std::string outName = contig + std::string(".fa");
    FILE* fp = std::fopen(outName.c_str(), "w");
    std::fwrite(">", sizeof(char), 1, fp);
    std::fwrite(contig, sizeof(char), std::strlen(contig), fp);
    std::fwrite("\n", sizeof(char), 1, fp);
    // split into 100nt/line
    int beg = 0;
    while(beg + 100 < len){
        fwrite(seq + beg, sizeof(char), 100, fp);
        std::fwrite("\n", sizeof(char), 1, fp);
        beg += 100;
    }
    if(beg < len - 1){
        fwrite(seq + beg, sizeof(char), len - beg, fp);
        std::fwrite("\n", sizeof(char), 1, fp);
    };
    std::fclose(fp);
    assert(fai_build(outName.c_str()) == 0);
    fai_destroy(fai);
}

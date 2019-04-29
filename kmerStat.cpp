#include "kseq.h"
#include <zlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

KSEQ_INIT(gzFile, gzread)

void statKmer(std::map<std::string, int32_t>& kcount, char* s, size_t l, int klen){
    for(size_t i = 0; i < l - klen + 1; ++i){
        std::string kmer(s + i, s + i + klen);
        if(kcount.find(kmer) != kcount.end()){
            ++kcount[kmer];
        }else{
            kcount[kmer] = 0;
        }
    }
}

int main(int argc, char** argv){
    if(argc <= 3){
        std::cout << argv[0] << " <infa> <kmerlen> <outf> " << std::endl;
        return 0;
    }
    char* infa = argv[1];
    int klen = std::atoi(argv[2]);
    char* outf = argv[3];
    std::ofstream fw(outf);
    std::map<std::string, int32_t> kcount;

    gzFile fp = gzopen(infa, "r");
    kseq_t* seq = kseq_init(fp);
    while(kseq_read(seq) >= 0){
        char* s = seq->seq.s;
        size_t l = seq->seq.l;
        statKmer(kcount, s, l, klen);
    }
    for(auto& e: kcount){
        fw << e.first << " " << e.second << "\n";
    }
    fw.close();
        
    gzclose(fp);
    kseq_destroy(seq);
}

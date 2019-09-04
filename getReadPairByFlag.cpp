#include "htslib/sam.h"
#include <kseq.h>
#include <zlib.h>
#include <string>
#include <iostream>
#include <unordered_set>
#include <iomanip>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv){
    if(argc < 7){
        std::cout << argv[0] << " <inbam> <infq1> <infq2> <flag> <outfq1> <outfq2> <drop>" << std::endl;
        return 0;
    }
    char* inbam = argv[1];
    char* infq1 = argv[2];
    char* infq2 = argv[3];
    uint16_t flag = std::atoi(argv[4]);
    char* outfq1 = argv[5];
    char* outfq2 = argv[6];
    bool drop = false;
    if(argc == 8){
        drop = true;
    }

    std::cout << "flag: " << flag << std::endl;
    std::cout << "drop: " << std::boolalpha << drop << std::endl;
    samFile* fp = sam_open(inbam, "r");
    bam_hdr_t* bh = sam_hdr_read(fp);
    bam1_t* br = bam_init1();
    std::unordered_set<std::string> rname;
    while(sam_read1(fp, bh, br) >= 0){
        if(br->core.flag & flag){
            if(!drop){
                rname.insert(bam_get_qname(br));
            }
        }else{
            if(drop){
                rname.insert(bam_get_qname(br));
            }
        }
    }
    std::cout << "should output read pairs: " << rname.size() << std::endl;
    bam_destroy1(br);
    bam_hdr_destroy(bh);
    sam_close(fp);

    gzFile fqi1 = gzopen(infq1, "r");
    gzFile fqi2 = gzopen(infq2, "r");
    kseq_t* seq1 = kseq_init(fqi1);
    kseq_t* seq2 = kseq_init(fqi2);
    gzFile fqo1 = gzopen(outfq1, "w");
    gzFile fqo2 = gzopen(outfq2, "w");
    gzsetparams(fqo1, 4, Z_DEFAULT_STRATEGY);
    gzsetparams(fqo2, 4, Z_DEFAULT_STRATEGY);
    int32_t outcount = 0;
    while(kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0){
        if(rname.count(seq1->name.s) > 0){
            gzputc(fqo1, '@');
            gzwrite(fqo1, seq1->name.s, seq1->name.l);
            gzputc(fqo1, '\n');
            gzwrite(fqo1, seq1->seq.s, seq1->seq.l);
            gzwrite(fqo1, "\n+\n", 3);
            gzwrite(fqo1, seq1->qual.s, seq1->qual.l);
            gzputc(fqo1, '\n');

            gzputc(fqo2, '@');
            gzwrite(fqo2, seq2->name.s, seq2->name.l);
            gzputc(fqo2, '\n');
            gzwrite(fqo2, seq2->seq.s, seq2->seq.l);
            gzwrite(fqo2, "\n+\n", 3);
            gzwrite(fqo2, seq2->qual.s, seq2->qual.l);
            gzputc(fqo2, '\n');
            ++outcount;
        }
    }
    gzflush(fqo1, Z_FINISH);
    gzflush(fqo2, Z_FINISH);
    gzclose(fqo1);
    gzclose(fqo2);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fqi1);
    gzclose(fqi2);
    fqo1 = NULL;
    fqo2 = NULL;
    fqi1 = NULL;
    fqi2 = NULL;
    std::cout << "total output read pairs: " << outcount << std::endl;
}

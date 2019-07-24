#include <htslib/faidx.h>
#include <iostream>
#include <string>
#include <util.h>

int main(int argc, char** argv){
    if(argc < 3){
        printf("%s <in.fa> <out.fa>\n", argv[0]);
        return 0;
    }
    faidx_t* fai = fai_load(argv[1]);
    int32_t seqlen = -1;
    char* chr16 = faidx_fetch_seq(fai, "chr16", 0, faidx_seq_len(fai, "chr16"), &seqlen);
    FILE* fp = fopen(argv[2], "w");
    fwrite(">simchr16\n", sizeof(char), 10, fp);
    // begin pos
    int beg = 1000000;
    fwrite(chr16 + beg - 3000, sizeof(char), 3000, fp); // get flank 3kb
    // inversion
    int inv1reg[2] = {beg, 1800 + beg};
    std::string inv1seq = std::string(chr16 + inv1reg[0], chr16 + inv1reg[1] + 1);
    util::reverseComplement(inv1seq);
    fwrite(inv1seq.c_str(), sizeof(char), inv1seq.length(), fp);// inv1reg 
    int inv2reg[2] = {6000 + beg, 8000 + beg};
    fwrite(chr16 + inv1reg[1] + 1, sizeof(char), inv2reg[0] - inv1reg[1] - 1, fp);// inter region
    std::string inv2seq = std::string(chr16 + inv2reg[0], chr16 + inv2reg[1] + 1);
    util::reverse(inv2seq);
    fwrite(inv2seq.c_str(), sizeof(char), inv2seq.length(), fp);// inv2reg
    // deletion
    int delreg1[2] = {9000 + beg, 10000 + beg};
    fwrite(chr16 + inv2reg[1] + 1, sizeof(char), delreg1[0] - inv2reg[1] - 1, fp);// inter region
    int delreg2[2] = {60000 + beg, 63000 + beg};
    fwrite(chr16 + delreg1[1] + 1, sizeof(char), delreg2[0] - delreg1[1] - 1, fp);// inter region
    int delreg3[2] = {69000 + beg, 70000 + beg};
    fwrite(chr16 + delreg2[1] + 1, sizeof(char), delreg3[0] - delreg2[1] - 1, fp);// inter region
    // duplication
    int dupreg1[2] = {80000 + beg, 82000 + beg};
    fwrite(chr16 + delreg3[1] + 1, sizeof(char), dupreg1[0] - delreg3[1] - 1, fp);// inter region
    fwrite(chr16 + dupreg1[0], sizeof(char), dupreg1[1] - dupreg1[0] + 1, fp);// dupreg1
    fwrite(chr16 + dupreg1[0], sizeof(char), dupreg1[1] - dupreg1[0] + 1, fp);// dupreg1
    int dupreg2[2] = {90000 + beg, 91000 + beg};
    fwrite(chr16 + dupreg1[1] + 1, sizeof(char), dupreg2[0] - dupreg1[1] - 1, fp);// inter region
    fwrite(chr16 + dupreg2[0], sizeof(char), dupreg2[1] - dupreg2[0] + 1, fp);// dupreg2
    fwrite(chr16 + dupreg2[0], sizeof(char), dupreg2[1] - dupreg2[0] + 1, fp);// dupreg2
    // translocation
    int trans1 = 100000 + beg;
    fwrite(chr16 + dupreg2[1] + 1, sizeof(char), trans1 - dupreg2[1], fp);// inter region
    free(chr16);
    char* chr17 = faidx_fetch_seq(fai, "chr17", 0, faidx_seq_len(fai, "chr17"), &seqlen);
    int trans2 = 200000 + beg;
    fwrite(chr17 + trans2, sizeof(char), 10000, fp);// trans region 10k flank
    free(chr17);
    // create faidx
    fclose(fp);
    assert(fai_build(argv[2]) == 0);
    fai_destroy(fai);
}

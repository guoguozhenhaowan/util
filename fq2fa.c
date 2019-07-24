#include <stdio.h>
#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv){
    if(argc < 3){
        printf("%s <in.fq> <out.fa>\n", argv[0]);
        return 0;
    }
    gzFile ifp = gzopen(argv[1], "r");
    kseq_t* seq = kseq_init(ifp);
    gzFile ofp = gzopen(argv[2], "w");
    while(kseq_read(seq) >= 0){
        gzwrite(ofp, ">", 1);
        gzwrite(ofp, seq->name.s, seq->name.l);
        gzwrite(ofp, "\n", 1);
        gzwrite(ofp, seq->seq.s, seq->seq.l);
        gzwrite(ofp, "\n", 1);
    }
    kseq_destroy(seq);
    gzclose(ifp);
    gzclose(ofp);
}

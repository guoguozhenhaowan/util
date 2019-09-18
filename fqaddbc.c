#include <stdio.h>
#include <kseq.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv){
    if(argc < 5){
        printf("%s <in.fq> <out.fq> <barcode> <quality>\n", argv[0]);
        return 0;
    }
    gzFile ifp = gzopen(argv[1], "r");
    kseq_t* seq = kseq_init(ifp);
    gzFile ofp = gzopen(argv[2], "w");
    char* bcseq = argv[3];
    int bclen = strlen(bcseq);
    char* quaseq = argv[4];
    int qualen = strlen(quaseq);
    while(kseq_read(seq) >= 0){
        gzputc(ofp, '@');
        gzwrite(ofp, seq->name.s, seq->name.l);
        gzsetparams(ofp, 1, Z_DEFAULT_STRATEGY);
        if(seq->comment.l){
            gzputc(ofp, ' ');
            gzwrite(ofp, seq->comment.s, seq->comment.l);
        }
        gzputc(ofp, '\n');
        gzwrite(ofp, bcseq, bclen);
        gzwrite(ofp, seq->seq.s, seq->seq.l);
        gzwrite(ofp, "\n+\n", 3);
        gzwrite(ofp, quaseq, qualen);
        gzwrite(ofp, seq->qual.s, seq->qual.l);
        gzputc(ofp, '\n');
    }
    kseq_destroy(seq);
    gzclose(ifp);
    gzclose(ofp);
}

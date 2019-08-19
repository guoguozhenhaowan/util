#include <kseq.h>
#include <zlib.h>
#include <stdio.h>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv){
    if(argc < 3){
        printf("%s <in.fq> <seq>\n", argv[0]);
        return 0;
    }
    char* infq = argv[1];
    char* qseq = argv[2];
    gzFile fp = gzopen(infq, "r");
    kseq_t* ks = kseq_init(fp);
    while(kseq_read(ks) >= 0){
        if(strstr(ks->seq.s, qseq)) printf("%s\n", ks->seq.s);
    }
    gzclose(fp);
    kseq_destroy(ks);
}

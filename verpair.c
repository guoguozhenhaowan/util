#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1;
    gzFile fp2;
	kseq_t *seq1;
    kseq_t *seq2;
	int l, m;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <r1.fasta> <r2.fasta> \n", argv[0]);
		return 1;
	}
	fp1 = gzopen(argv[1], "r");
    fp2 = gzopen(argv[2], "r");
	seq1 = kseq_init(fp1);
    seq2 = kseq_init(fp2);
	while ((l = kseq_read(seq1)) >= 0 && (m = kseq_read(seq2)) >= 0) {
        if(strcmp(seq1->name.s, seq2->name.s)){
            printf("name: %s\t%s\n", seq1->name.s, seq2->name.s);
        }
    }
	kseq_destroy(seq1);
    kseq_destroy(seq2);
	gzclose(fp1);
    gzclose(fp2);
	return 0;
}
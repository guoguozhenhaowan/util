#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <kseq.h>
#include <util.h>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1;
    gzFile fp2;
	kseq_t *seq1;
    kseq_t *seq2;
	int l, m;
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <fqr1> <fqr2> <readname>\n", argv[0]);
		return 1;
	}
	fp1 = gzopen(argv[1], "r");
    fp2 = gzopen(argv[2], "r");
	seq1 = kseq_init(fp1);
    seq2 = kseq_init(fp2);
    std::cout << "r1beg: " << util::currentTime() << std::endl;
    while((l = kseq_read(seq1)) >= 0);
    std::cout << "r1end: " << util::currentTime() << std::endl;
    std::cout << "r2beg: " << util::currentTime() << std::endl;
    while((m = kseq_read(seq2)) >= 0);
    std::cout << "r2end: " << util::currentTime() << std::endl;
	kseq_destroy(seq1);
    kseq_destroy(seq2);
	gzclose(fp1);
    gzclose(fp2);
	return 0;
}

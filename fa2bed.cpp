#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <kseq.h>
#include <CLI.hpp>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        std::cout << "fa2bed -h for help" << std::endl;
        return 0;
    }
    std::string infa;
    std::string outbed;
    
    CLI::App app{"fasta to bed"};
    app.add_option("-i,--in", infa, "input fasta file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-o,--out", outbed, "output bed file");
    CLI_PARSE(app, argc, argv);

    std::string bedrec = "";
    gzFile fpr = gzopen(infa.c_str(), "r");
    kseq_t* seq = kseq_init(fpr);
    int l;
    while((l = kseq_read(seq)) >= 0)
    {
        bedrec.append(seq->name.s);
        bedrec.append("\t0\t");
        bedrec.append(std::to_string(seq->seq.l));
        bedrec.append("\n");
    }

    std::ofstream fw;
    fw.open(outbed.c_str(), std::ios::out);
    fw << bedrec;
    fw.close();
    kseq_destroy(seq);
    return 0;
}

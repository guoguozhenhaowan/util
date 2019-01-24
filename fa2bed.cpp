#include <iostream>
#include <fstream>
#include <string>
#include <libgen.h>
#include <zlib.h>
#include <kseq.h>
#include <CLI.hpp>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv)
{
    std::string sys_cmd = std::string(argv[0]) + " -h";
    if(argc < 2){std::system(sys_cmd.c_str()); return 0;}
    
    std::string version = "0.0.0";
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    std::string infa;
    std::string outbed;
    
    CLI::App app{"program: " + std::string(basename(argv[0])) + "\nversion: " + version + "\nupdated: " + cmp_time};
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

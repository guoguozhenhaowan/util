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
    std::string outfa;
    
    CLI::App app{"program: " + std::string(basename(argv[0])) + "\nversion: " + version + "\nupdated: " + cmp_time};
    app.add_option("-i,--in", infa, "input fasta file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-o,--out", outfa, "output fasta file");
    CLI_PARSE(app, argc, argv);

    std::string farec = "";
    gzFile fpr = gzopen(infa.c_str(), "r");
    kseq_t* seq = kseq_init(fpr);
    std::ofstream fw(outfa);
    int l;
    while((l = kseq_read(seq)) >= 0)
    {
        farec.clear();
        farec.append(">");
        farec.append(seq->name.s);
        farec.append("\n");
        size_t index = 0;
        while(index + 50 <= seq->seq.l){
            farec.append(std::string(seq->seq.s + index, 50));
            farec.append("\n");
            index += 50;
        }
        if(index < seq->seq.l){
            farec.append(seq->seq.s + index);
            farec.append("\n");
        }
        farec.append("\n");
        fw << farec;
    }
    fw.close();
    kseq_destroy(seq);
    gzclose(fpr);
    return 0;
}

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
        std::cout << argv[0] << " -h for help" << std::endl;
        return 0;
    }
    std::string infa;
    std::string patt;
    std::string outfa = "ext.fa";
    unsigned long long maxline = 1000000ULL;
    
    CLI::App app{"fasta extraction"};
    app.add_option("-i,--in", infa, "input fasta file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-p,--patt", patt, "fixed pattern to search")->required(true);
    CLI::Option* pout = app.add_option("-o,--out", outfa, "output fasta file");
    app.add_option("-m,--max", maxline, "maxline in memory", true);
    CLI_PARSE(app, argc, argv);

    gzFile fpr = gzopen(infa.c_str(), "r");
    kseq_t* seq = kseq_init(fpr);
    std::ofstream fw;
    fw.open(outfa.c_str(), std::ios::out | std::ios::binary);
    std::string tmp_str = "";
    std::string tmp_nam;
    std::string tmp_com;
    tmp_str.reserve(maxline * 2000);
    int l = 1;
    while(l >= 0)
    {
        tmp_str.clear();
        for(unsigned long long i = 0; i < maxline; ++i)
        {
            l = kseq_read(seq);
            if(l >= 0)
            {
                tmp_nam = seq->name.s;
                tmp_com = seq->comment.s;
                if((tmp_nam.find(patt) != std::string::npos) || (tmp_com.find(patt) != std::string::npos))
                {
                    tmp_str.append(">");
                    tmp_str.append(tmp_nam);
                    tmp_str.append("\n");
                    tmp_str.append(seq->seq.s);
                    tmp_str.append("\n");
                }
            }
            else
            {
                break;
            }
        }
        fw << tmp_str;
    }
    fw.close();
    kseq_destroy(seq);
    return 0;
}

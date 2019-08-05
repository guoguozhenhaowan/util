#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <stdio.h>
#include <util.h>
#include <map>
#include <fstream>

int main(int argc, char** argv){
    if(argc < 4){
        printf("%s <refGene.sorted.gtf.gz> <trans2gene.tsv> <refgene.modified.sorted.gtf.gz>\n", argv[0]);
        return 0;
    }
    // parse transcript to genename relations
    std::map<std::string, std::string> tr2gn;
    std::vector<std::string> vstr;
    std::ifstream fr(argv[2]);
    std::string tmpStr;
    while(std::getline(fr, tmpStr)){
        util::split(tmpStr, vstr, "\t");
        tr2gn[vstr[0]] = vstr[1];
    }
    fr.close();
    // modify gtf group information
    BGZF* ofp = bgzf_open(argv[3], "wb");
    BGZF* ifp = bgzf_open(argv[1], "rb");
    kstring_t str = {0, 0, 0};
    std::vector<std::string> gvstr;
    while(bgzf_getline(ifp, '\n', &str) >= 0){
        util::split(str.s, vstr, "\t");
        util::split(vstr[8], gvstr, " ");
        std::string trsid = util::replace(gvstr[1], "\"", "");
        trsid = util::replace(trsid, ";", "");
        gvstr[0] = "gene_name";
        gvstr[1] = "\"" + tr2gn[trsid] + "\";";
        std::string gnstr;
        util::join(gvstr, gnstr, " ");
        vstr[8] = gnstr;
        std::string rec;
        util::join(vstr, rec, "\t");
        rec.append("\n");
        assert(bgzf_write(ofp, rec.c_str(), rec.size()) >= 0);
    }
    bgzf_close(ofp);
    bgzf_close(ifp);
    // create index
    tbx_conf_t tc = tbx_conf_gff;
    tbx_index_build(argv[3], 0, &tc);
    // test output file
    htsFile* tfp = hts_open(argv[3], "r");
    tbx_t* tbx = tbx_index_load(argv[3]);
    hts_itr_t* itr = tbx_itr_queryi(tbx, tbx_name2id(tbx, "chr1"), 34611, 36800);
    while(tbx_itr_next(tfp, tbx, itr, &str) >= 0){
        std::cout << str.s << std::endl;
    }
    hts_close(tfp);
    tbx_itr_destroy(itr);
    tbx_destroy(tbx);
}

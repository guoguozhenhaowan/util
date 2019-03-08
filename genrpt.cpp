/* generate report file
 * required files listed below
 * lib.filt.log         : rna filter log file
 * lib_DDP_QC.tsv       : DDP QC tsv file
 * lib_IDP_QC.tsv       : IDP QC tsv file
 * lib.FusionReport.txt : fusionMap result
 * lib/abundance.tsv    : kallisto result
 * ensebml2genename     : NCBI ensembl to genename file
 * geneset.tsv          : geneset tsv file
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <utility>
#include <tuple>
#include <map>
#include "CLI.hpp"
#include "xlsxwriter.h"

namespace genrpt{
    struct args{
        std::string filtlog; // lib.filt.log
        std::string ddpqc; // lib_DDP_QC.tsv
        std::string idpqc; // lib_IDP_QC.tsv
        std::string fusionrpt; // lib.FusionReport.txt
        std::string kabundance; // lib/abundance.tsv
        std::string ens2gname; // ensebml2genename
        std::string gset; // geneset.tsv
        std::string outf; // lib.report
    };
    
    struct stats{
        std::map<std::string, double> gexp; // gene expression of gset 
        std::map<std::string, double> gout; // other gene expression 
        double rratio = 0.0; // rrna ratio 
        double rfusion = 0.0; // reads support fusion
    };
    
    void get_exp(genrpt::stats& s, const genrpt::args& a);
    void get_fratio(genrpt::stats& s, const genrpt::args& a);
    double get_rratio(const std::string& filtlog);
    void gen_finarpt(const genrpt::stats& s, const genrpt::args& a);
    void split_line(std::vector<std::string>& v, const std::string& str, const char& sep){
        v.clear();
        std::stringstream ss(str);
        std::string tmpstr;
        while(std::getline(ss, tmpstr, sep)){
            v.push_back(tmpstr);
        }
    }
}

int main(int argc, char** argv){
    genrpt::args a;
    CLI::App app("generate report program");
    app.add_option("-f", a.filtlog, "lib.filt.log")->required()->check(CLI::ExistingFile);
    app.add_option("-d", a.ddpqc, "lib_DDP_QC.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-i", a.idpqc, "lib_IDP_QC.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-s", a.fusionrpt, "lib.FusionReport.txt")->required()->check(CLI::ExistingFile);
    app.add_option("-k", a.kabundance, "lib/abundance.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-e", a.ens2gname, "ensebml2genename")->required()->check(CLI::ExistingFile);
    app.add_option("-g", a.gset, "geneset.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-o", a.outf, "lib.report")->required();
    CLI_PARSE(app, argc, argv);

    genrpt::stats s;
    s.rratio = genrpt::get_rratio(a.filtlog);
    genrpt::get_exp(s, a);
    genrpt::get_fratio(s, a);
    genrpt::gen_finarpt(s, a);
}

void genrpt::get_exp(genrpt::stats& s, const genrpt::args& a){
    std::map<std::string, std::vector<std::string>> gem; // gene to ensembl map
    std::map<std::string, double> eex; // ensembl express 
    std::ifstream fe2g(a.ens2gname), fexp(a.kabundance), fgset(a.gset);
    std::stringstream ss1, ss2, ss3;
    ss1 << fe2g.rdbuf();
    ss2 << fexp.rdbuf();
    ss3 << fgset.rdbuf();
    std::string ens, gene, tid, line; 
    double l, el, ec, tpm;
    std::getline(ss1, line);
    std::getline(ss2, line);
    while(ss1 >> ens >> gene){
        if(gem.find(gene) == gem.end()){
            gem[gene] = {ens};
        }else{
            gem[gene].push_back(ens);
        }
    }
    while(ss2 >> tid >> l >> el >> ec >> tpm){
        eex[tid.substr(0, tid.find_first_of("."))] = tpm;
    }
    for(auto& e: gem){
        double exp = 0.0;
        for(auto& f: e.second){
            exp += ((eex.find(f) != eex.end()) ? eex[f] : 0);
        }
        if(exp > 0){
            s.gout[e.first] = exp;
        }
    }
    std::vector<std::string> giset;
    std::vector<std::string> vsg;
    while(std::getline(ss3, line)){
        genrpt::split_line(vsg, line, '\t');
        gene = vsg[0];
        if(s.gout.find(gene) != s.gout.end()){
            s.gexp[gene] = s.gout[gene];
        }else{
            s.gexp[gene] = 0;
        }
        giset.push_back(gene);
    }
    for(auto& e: giset){
        s.gout.erase(e);
    }
}

void genrpt::get_fratio(genrpt::stats& s, const genrpt::args& a){
    std::pair<std::string, std::string> tpair;
    std::vector<std::string> tfus(26);
    std::vector<std::string> vtfus; 
    std::ifstream frq(a.ddpqc), frf(a.fusionrpt);
    std::stringstream ss1, ss2;
    ss1 << frq.rdbuf();
    ss2 << frf.rdbuf();
    std::string line;
    size_t ttreads;
    std::vector<std::string> vs;
    while(std::getline(ss1, line)){
        genrpt::split_line(vs, line, '\t');
        tpair.first = vs[0];
        tpair.second = vs[1];
        if(tpair.first == "Total Reads"){
            ttreads = std::atoi(tpair.second.c_str());
            break;
        }
    }
    std::getline(ss2, line);
    std::vector<size_t> vpos = {0};
    while(std::getline(ss2, line)){
        genrpt::split_line(vs, line, '\t');
        std::string seedc = vs[2];
        std::string resuc = vs[3];
        s.rfusion += std::atof(seedc.c_str()) + std::atof(resuc.c_str());
    }
    s.rfusion /= ttreads;
}

double genrpt::get_rratio(const std::string& filtlog){
    std::ifstream fr(filtlog);
    std::string line, suffix;
    std::vector<std::string> vs;
    std::stringstream ss;
    while(std::getline(fr, line)){
        if(line[0] == '%'){
            genrpt::split_line(vs, line, ' ');
            suffix = vs[vs.size() - 1];
            return(std::atof(suffix.c_str()));
        }
    }
    return 0;
}

void genrpt::gen_finarpt(const genrpt::stats& s, const genrpt::args& a){
    lxw_workbook* workbook = new_workbook(a.outf.c_str());
    // sheet IDP
    lxw_worksheet* sheet = workbook_add_worksheet(workbook, "IDP");
    std::pair<std::string, std::string> dpair; 
    std::string line, prefix, suffix;
    std::ifstream fr(a.idpqc);
    int row = 0, col = 0;
    size_t maxc1len = 0, maxc2len = 0;
    while(std::getline(fr, line)){
        auto p = line.find("\t");
        prefix = line.substr(0, p);
        suffix = line.substr(p + 1);
        worksheet_write_string(sheet, row, 0, prefix.c_str(), NULL);
        worksheet_write_string(sheet, row++, 1, suffix.c_str(), NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    fr.close();
    // sheet DDP
    sheet = workbook_add_worksheet(workbook, "DDP");
    fr.open(a.ddpqc);
    row = 0, col = 0;
    maxc1len = 0, maxc2len = 0;
    while(std::getline(fr, line)){
        auto p = line.find("\t");
        prefix = line.substr(0, p);
        suffix = line.substr(p + 1);
        worksheet_write_string(sheet, row, 0, prefix.c_str(), NULL);
        worksheet_write_string(sheet, row++, 1, suffix.c_str(), NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    fr.close();
    // sheet GSETEXP
    sheet = workbook_add_worksheet(workbook, "GSETEXP");
    row = 0, col = 0;
    maxc1len = 0, maxc2len = 0;
    worksheet_write_string(sheet, row, 0, "Gene", NULL);
    worksheet_write_string(sheet, row++, 1, "TPM", NULL);
    std::ostringstream oss;
    for(auto& e: s.gexp){
        prefix = e.first;
        oss.clear();
        oss.str("");
        oss << e.second;
        suffix = oss.str();
        worksheet_write_string(sheet, row, 0, prefix.c_str(), NULL);
        worksheet_write_string(sheet, row++, 1, suffix.c_str(), NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    // sheet GOSETEXP
    sheet = workbook_add_worksheet(workbook, "GOSETEXP");
    row = 0, col = 0;
    maxc1len = 0, maxc2len = 0;
    worksheet_write_string(sheet, row, 0, "Gene", NULL);
    worksheet_write_string(sheet, row++, 1, "TPM", NULL);
    for(auto& e: s.gout){
        prefix = e.first;
        oss.clear();
        oss.str("");
        oss << e.second;
        suffix = oss.str();
        worksheet_write_string(sheet, row, 0, prefix.c_str(), NULL);
        worksheet_write_string(sheet, row++, 1, suffix.c_str(), NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    // sheet fusionmap
    sheet = workbook_add_worksheet(workbook, "FusionMap");
    row = 0, col = 0;
    maxc1len = 0, maxc2len = 0;
    std::vector<size_t> vclen;
    std::vector<std::string> vrec;
    fr.open(a.fusionrpt.c_str());
    std::getline(fr, line);
    genrpt::split_line(vrec, line, '\t');
    for(size_t i = 0; i < vrec.size(); ++i){
        vclen.push_back(vrec[i].length());
        worksheet_write_string(sheet, row, col++, vrec[i].c_str(), NULL);
    }
    while(std::getline(fr, line)){
        ++row;
        col = 0;
        genrpt::split_line(vrec, line, '\t');
        for(size_t i = 0; i < vrec.size(); ++i){
            vclen[i] = std::max(vclen[i], vrec[i].length());
            worksheet_write_string(sheet, row, col++, vrec[i].c_str(), NULL);
        }
    }
    for(size_t i = 0; i < vclen.size(); ++i){
        worksheet_set_column(sheet, i, i, vclen[i], NULL);
    }
    // sheet Extra
    row = 0, col = 0;
    maxc1len = 0,  maxc2len = 0;
    sheet = workbook_add_worksheet(workbook, "Extra");
    // rrna ratio
    oss.clear();
    oss.str("");
    oss << s.rratio;
    suffix = oss.str();
    maxc2len = std::max(suffix.length(), maxc2len);
    worksheet_write_string(sheet, row, col++, "rRNA Ratio", NULL);
    worksheet_write_string(sheet, row++, col++, suffix.c_str(), NULL);
    // fusion reads ratio
    col = 0;
    oss.clear();
    oss.str("");
    oss << s.rfusion;
    suffix = oss.str();
    maxc2len = std::max(suffix.length(), maxc2len);
    col = 0; 
    worksheet_write_string(sheet, row, col++, "FusionReadsRatio", NULL);
    worksheet_write_string(sheet, row++, col++, suffix.c_str(), NULL);
    worksheet_set_column(sheet, 0, 0, 16, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    workbook_close(workbook);
}

/* generate report file
 * required files listed below
 * lib.fqtool.json      : fqtool json statistic file
 * lib.filter.json      : filter json statistic file
 * lib.BasicQC.json     : bamqc json statistic file
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
#include <cmath>
#include <utility>
#include <algorithm>
#include <tuple>
#include <map>
#include <CLI.hpp>
#include <json.hpp>
#include "xlsxwriter.h"

namespace genrpt{
    struct args{
        std::string fqqc;       ///< lib.fqtool.json
        std::string filtlog;    ///< lib.filter.json
        std::string bamqc;      ///< lib.BasicQC.json
        std::string fusionrpt;  ///< lib.FusionReport.txt
        std::string kabundance; ///< lib/abundance.tsv
        std::string ens2gname;  ///< ensebml2genename
        std::string gset;       ///< geneset.tsv
        std::string outf;       ///< lib.report
    };
    
    struct stats{
        std::map<std::string, double> gexp; ///< gene expression of gset
        std::map<std::string, double> gout; ///< other gene expression
        double rratio = 0.0;                ///< rrna ratio
        double rfusion = 0.0;               ///< reads support fusion
    };
    
    void get_exp(genrpt::stats& s, const genrpt::args& a);
    void get_fratio(genrpt::stats& s, const genrpt::args& a);
    void gen_fqqcsheet(lxw_worksheet* sheet, jsn::json& j);
    void gen_ncrnasheet(lxw_worksheet* sheet, jsn::json& j);
    void gen_bamqcsheet(lxw_worksheet* sheet, jsn::json& j);
    void gen_finarpt(const genrpt::stats& s, const genrpt::args& a);
    inline bool starts_with(const std::string& str, const std::string& pre){
        if(str.length() < pre.length()){
            return false;
        }else{
            return std::equal(pre.begin(), pre.end(), str.begin());
        }
    }
    inline void split_line(std::vector<std::string>& vec, const std::string& str, std::string sep = " "){
        vec.clear();
        std::string::size_type las, cur;
        las = cur = str.find_first_not_of(sep);
        while((las = str.find_first_not_of(sep, las)) != std::string::npos){
            cur = str.find_first_of(sep, las);
            if(cur != std::string::npos){
                vec.push_back(str.substr(las, cur - las));
            }else{
                vec.push_back(str.substr(las));
                break;
            }
            las = cur;
        }
    }
    inline std::string replace(const std::string& str, const std::string& pat, const std::string& des){
        std::string ret;
        std::string::size_type las = 0, cur = 0;
        while((cur = str.find(pat, cur)) != std::string::npos){
            ret.append(str.substr(las, cur - las));
            ret.append(des);
            cur += pat.length();
            las = cur;
        }
        if(las != std::string::npos){
            ret.append(str.substr(las));
        }
        return ret;
    }

}

int main(int argc, char** argv){
    genrpt::args a;
    CLI::App app("generate report program");
    app.add_option("-q", a.fqqc, "lib.fqtool.json")->required()->check(CLI::ExistingFile);
    app.add_option("-f", a.filtlog, "lib.filter.json")->required()->check(CLI::ExistingFile);
    app.add_option("-b", a.bamqc, "lib.bamqc.json")->required()->check(CLI::ExistingFile);
    app.add_option("-s", a.fusionrpt, "lib.FusionReport.txt")->required()->check(CLI::ExistingFile);
    app.add_option("-k", a.kabundance, "lib/abundance.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-e", a.ens2gname, "ensebml2genename")->required()->check(CLI::ExistingFile);
    app.add_option("-g", a.gset, "geneset.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-o", a.outf, "lib.report")->required();
    CLI_PARSE(app, argc, argv);

    genrpt::stats s;
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
        genrpt::split_line(vsg, line, "\t");
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
    // get total reads
    std::ifstream frq(a.bamqc);
    jsn::json jqc;
    frq >> jqc;
    size_t ttreads = jqc["DupIncludedQC"]["TotalReads"];
    frq.close();
    // get fusion reads
    std::ifstream frf(a.fusionrpt);
    std::vector<std::string> vs;
    std::stringstream ss;
    ss << frf.rdbuf();
    std::string line;
    std::getline(ss, line);
    while(std::getline(ss, line)){
        genrpt::split_line(vs, line, "\t");
        std::string seedc = vs[2];
        std::string resuc = vs[3];
        s.rfusion += std::atof(seedc.c_str()) + std::atof(resuc.c_str());
    }
    // get fusion reads ratio
    s.rfusion /= ttreads;
}

void genrpt::gen_fqqcsheet(lxw_worksheet* sheet, jsn::json& j){
    size_t row = 0, maxc1len = 0, maxc2len = 0;
    std::vector<std::string> items = {"TotalReads", "TotalBases", "Read1Length", "Read2Length", "GCRate", "Q20Bases", "Q20BaseRate", "Q30Bases", "Q30BaseRate"};
    std::stringstream ss;
    std::string colStr;
    for(uint32_t i = 0; i < items.size(); ++i){
        ss.clear();
        ss.str("");
        ss << j[items[i]];
        colStr = ss.str();
        colStr = replace(colStr, "\"", "");
        worksheet_write_string(sheet, row, 0, items[i].c_str(), NULL);
        worksheet_write_number(sheet, row++, 1, j[items[i]], NULL);
        maxc1len = std::max(maxc1len, items[i].length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
}

void genrpt::gen_ncrnasheet(lxw_worksheet* sheet, jsn::json& j){
    size_t row = 0, maxc1len = 0, maxc2len = 0;
    std::vector<std::string> genItems = {"ReadsIn", "ReadsGot", "ReadsDrop", "DropRate"};
    std::stringstream ss;
    std::string colStr;
    jsn::json jgen = j["Summary"];
    for(uint32_t i = 0; i < genItems.size(); ++i){
        ss.clear();
        ss.str("");
        ss << jgen[genItems[i]];
        colStr = ss.str();
        colStr = replace(colStr, "\"", "");
        worksheet_write_string(sheet, row, 0, genItems[i].c_str(), NULL);
        worksheet_write_number(sheet, row++, 1, jgen[genItems[i]], NULL);
        maxc1len = std::max(maxc1len, genItems[i].length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    std::unordered_map<std::string, int32_t> jmap = j["FilterCount"];
    for(auto& e: jmap){
        ss.clear();
        ss.str("");
        ss << e.second;
        colStr = ss.str();
        colStr = replace(colStr, "\"", "");
        worksheet_write_string(sheet, row, 0, e.first.c_str(), NULL);
        worksheet_write_number(sheet, row++, 1, e.second, NULL);
        maxc1len = std::max(maxc1len, e.first.length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);

}
void genrpt::gen_bamqcsheet(lxw_worksheet* sheet, jsn::json& j){
    size_t row = 0, maxc1len = 0, maxc2len = 0;
    std::vector<std::string> genItems = {"ReadLength", "TotalReads", "RegionSize", "InsertSize", "EffectiveReads", "EffectiveReadsRate", 
                                      "EffectiveBases", "EffectiveBasesRate", "MismatchBases", "MismatchBasesRate", "UniqMappedReads", 
                                      "UniqMappedReadsRate", "PassedLowQCReads", "PassedLowQCReadsRate", "OnTargetReads", "OnTargetReadsRate",
                                      "OnTargetBases", "OnTargetBasesRate", "OnTargetMismatches", "OnTargetMismatchesRate"};
    std::vector<std::string> covItems = {"1XTargetRegionCoverage", "20XTargetRegionCoverage", "30XTargetRegionCoverage", "50XTargetRegionCoverage", 
                                      "100XTargetRegionCoverage", "200XTargetRegionCoverage", "300XTargetRegionCoverage", "500XTargetRegionCoverage", 
                                      "1000XTargetRegionCoverage", "2000XTargetRegionCoverage", "3000XTargetRegionCoverage", "5000XTargetRegionCoverage",
                                      "8000XTargetRegionCoverage", "10000XTargetRegionCoverage", "20000XTargetRegionCoverage", "30000XTargetRegionCoverage", 
                                      "50000XTargetRegionCoverage"};
    std::stringstream ss;
    std::string colStr;
    for(uint32_t i = 0; i < genItems.size(); ++i){
        ss.clear();
        ss.str("");
        ss << j[genItems[i]];
        colStr = ss.str();
        colStr = replace(colStr, "\"", "");
        worksheet_write_string(sheet, row, 0, genItems[i].c_str(), NULL);
        if(genItems[i] == "InsertSize"){
            worksheet_write_string(sheet, row++, 1, colStr.c_str(), NULL);
        }else{
            worksheet_write_number(sheet, row++, 1, j[genItems[i]], NULL);
        }
        maxc1len = std::max(maxc1len, genItems[i].length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    jsn::json jcov = j["Coverage"];
    for(uint32_t i = 0; i < covItems.size(); ++i){
        ss.clear();
        ss.str("");
        ss << jcov[covItems[i]];
        colStr = ss.str();
        worksheet_write_string(sheet, row, 0, covItems[i].c_str(), NULL);
        worksheet_write_number(sheet, row++, 1, jcov[covItems[i]], NULL);
        maxc1len = std::max(maxc1len, covItems[i].length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
}

void genrpt::gen_finarpt(const genrpt::stats& s, const genrpt::args& a){
    lxw_workbook* workbook = new_workbook(a.outf.c_str());
    std::ifstream fr;
    // sheet FastqQC
    jsn::json jfqqc;
    fr.open(a.fqqc.c_str());
    fr >> jfqqc;
    lxw_worksheet* sheet = workbook_add_worksheet(workbook, "FastqQC");
    gen_fqqcsheet(sheet, jfqqc["Summary"]["AfterFiltering"]);
    fr.close();
    // bamqc sheets
    jsn::json jbamqc;
    fr.open(a.bamqc.c_str());
    fr >> jbamqc;
    // sheet DupIncludeQC
    sheet = workbook_add_worksheet(workbook, "DupIncludeQC");
    gen_bamqcsheet(sheet, jbamqc["DupIncludedQC"]);
    // sheet DupExcludeQC
    sheet = workbook_add_worksheet(workbook, "DupExcludeQC");
    gen_bamqcsheet(sheet, jbamqc["DupExcludeQC"]);
    fr.close();
    // sheet OnTargetGeneExp
    sheet = workbook_add_worksheet(workbook, "OnTargetGeneExp");
    int row = 0, col = 0;
    size_t maxc1len = 0, maxc2len = 0;
    std::string prefix, suffix, line;
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
        worksheet_write_number(sheet, row++, 1, e.second, NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    // sheet AllGeneExp
    sheet = workbook_add_worksheet(workbook, "AllGeneExp");
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
        worksheet_write_number(sheet, row++, 1, e.second, NULL);
        maxc1len = std::max(maxc1len, prefix.length());
        maxc2len = std::max(maxc2len, suffix.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    // sheet FusionResult
    sheet = workbook_add_worksheet(workbook, "FusionResult");
    row = 0, col = 0;
    maxc1len = 0, maxc2len = 0;
    std::vector<size_t> vclen;
    std::vector<std::string> vrec;
    fr.open(a.fusionrpt.c_str());
    std::getline(fr, line);
    genrpt::split_line(vrec, line, "\t");
    for(size_t i = 0; i < vrec.size(); ++i){
        vclen.push_back(vrec[i].length());
        worksheet_write_string(sheet, row, col++, vrec[i].c_str(), NULL);
    }
    while(std::getline(fr, line)){
        ++row;
        col = 0;
        genrpt::split_line(vrec, line, "\t");
        for(size_t i = 0; i < vrec.size(); ++i){
            vclen[i] = std::max(vclen[i], vrec[i].length());
            worksheet_write_string(sheet, row, col++, vrec[i].c_str(), NULL);
        }
    }
    for(size_t i = 0; i < vclen.size(); ++i){
        worksheet_set_column(sheet, i, i, vclen[i], NULL);
    }
    fr.close();
    // sheet NCRnaRatio
    fr.open(a.filtlog.c_str());
    jsn::json jfilter;
    fr >> jfilter;
    sheet = workbook_add_worksheet(workbook, "NCRnaRatio");
    gen_ncrnasheet(sheet, jfilter["FilterResult"]);
    fr.close();
    // sheet Extra
    row = 0, col = 0;
    maxc1len = 0,  maxc2len = 0;
    sheet = workbook_add_worksheet(workbook, "Extra");
    oss.clear();
    oss.str("");
    oss << s.rfusion;
    suffix = oss.str();
    maxc2len = std::max(suffix.length(), maxc2len);
    worksheet_write_string(sheet, row, col++, "FusionReadsRatio", NULL);
    worksheet_write_number(sheet, row, col++, s.rfusion, NULL);
    worksheet_set_column(sheet, 0, 0, 16, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    workbook_close(workbook);
}

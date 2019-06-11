/* generate report file
 * required files listed below
 * lib.filter.json      : filter json statistic file
 * lib.BasicQC.json     : bamqc json statistic file
 * lib/abundance.tsv    : kallisto result
 * ensebml2genename     : NCBI ensembl to genename file
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
#include "util.h"
#include "xlsxwriter.h"

namespace genrpt{
    /** class to store arguments input */
    struct Options{
        std::string filterResult;       ///< lib.filter.json
        std::string bamqcResult;        ///< lib.BasicQC.json
        std::string kallistoResult;     ///< lib/abundance.tsv
        std::string ensemblToGeneTable; ///< ensebml2genename
        std::string outReportFile;      ///< lib.report
    };
    
    /** class to store various statistical items */
    struct Stats{
        std::map<std::string, double> geneExpress; ///< gene expression from kallisto
    };
   
    /** generate gene expression level sheet
     * @param s reference of Stats object
     * @param o reference of Options object
     */
    void getGeneExpress(genrpt::Stats& s, const genrpt::Options& o);
    
    /** generate rna classification sheet
     * @param sheet pointer to rna classification lxw_worksheet
     * @param jfilter reference of filter result 
     * @param jbamqc reference of bamqc result
     */
    void genRNAClassSheet(lxw_worksheet* sheet, jsn::json& jfilter, jsn::json& jbamqc);
    
    /** generate bamqc sheet
     * @param sheet pointer to bamqc lxw_worksheet
     * @param j pointer to json object storing bamqc result
     */
    void genBamQCSheet(lxw_worksheet* sheet, jsn::json& j);

    /** generate 3'utr length
     * @param sheet pointer to 3utrlen lxw_worksheet
     * @param jbamqc pointer to json object storing bamqc result
     */
    void gen3UTRLenSheet(lxw_worksheet* sheet, jsn::json& jbamqc);
    
    /** generate final report 
     * @param s reference of Stats object
     * @param a reference of Options object
     */
    void genFinalReport(const genrpt::Stats& s, const genrpt::Options& a);
}

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    genrpt::Options opt;
    CLI::App app("generate report program");
    app.add_option("-f", opt.filterResult, "lib.filter.json")->required()->check(CLI::ExistingFile);
    app.add_option("-b", opt.bamqcResult, "lib.bamqc.json")->required()->check(CLI::ExistingFile);
    app.add_option("-k", opt.kallistoResult, "lib/abundance.tsv")->required()->check(CLI::ExistingFile);
    app.add_option("-e", opt.ensemblToGeneTable, "ensebml2genename")->required()->check(CLI::ExistingFile);
    app.add_option("-o", opt.outReportFile, "lib.report")->required();
    CLI_PARSE(app, argc, argv);

    genrpt::Stats s;
    genrpt::getGeneExpress(s, opt);
    genrpt::genFinalReport(s, opt);
}

void genrpt::getGeneExpress(genrpt::Stats& s, const genrpt::Options& opt){
    std::map<std::string, std::vector<std::string>> gem; // gene to ensembl map
    std::map<std::string, double> eex; // ensembl express 
    std::ifstream fe2g(opt.ensemblToGeneTable), fexp(opt.kallistoResult);
    std::stringstream ss1, ss2;
    ss1 << fe2g.rdbuf();
    ss2 << fexp.rdbuf();
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
            s.geneExpress[e.first] = exp;
        }
    }
}

void genrpt::genRNAClassSheet(lxw_worksheet* sheet, jsn::json& jfilter, jsn::json& jbamqc){
    // stat items
    std::vector<std::string> statKeys = {"TotalReads", "mRNA", "lincRNA", "snRNA", "snoRNA", "rRNA", "miRNA", "Others"};
    std::map<std::string, int32_t> statMap;
    int32_t classifiedReads = 0;
    jsn::json jCount = jfilter["FilterCount"];
    statMap["TotalReads"] = jfilter["Summary"]["ReadsIn"];
    statMap["mRNA"] = jbamqc["EffectiveReads"];
    statMap["mRNA"] /= 2;
    classifiedReads += statMap["mRNA"];
    //lincRNA
    if(jCount.find("lincRNA") == jCount.end()){
        statMap["lincRNA"] = 0;
    }else{
        statMap["lincRNA"] = jCount["lincRNA"];
    }
    classifiedReads += statMap["lincRNA"];
    //snRNA
    if(jCount.find("snRNA") == jCount.end()){
        statMap["snRNA"] = 0;
    }else{
        statMap["snRNA"] = jCount["snRNA"];
    }
    classifiedReads += statMap["snRNA"];
    //snoRNA
    if(jCount.find("snoRNA") == jCount.end()){
        statMap["snoRNA"] = 0;
    }else{
        statMap["snoRNA"] = jCount["snoRNA"];
    }
    classifiedReads += statMap["snoRNA"];
    //rRNA
    if(jCount.find("rRNA") == jCount.end()){
        statMap["rRNA"] = 0;
    }else{
        statMap["rRNA"] = jCount["rRNA"];
    }
    classifiedReads += statMap["rRNA"];
    //miRNA
    if(jCount.find("miRNA") == jCount.end()){
        statMap["miRNA"] = 0;
    }else{
        statMap["miRNA"] = jCount["miRNA"];
    }
    classifiedReads += statMap["miRNA"];
    //Others
    statMap["Others"] = statMap["TotalReads"] - classifiedReads;
    //output
    std::stringstream ss;
    std::string colStr;
    size_t row = 0, maxc1len = 0, maxc2len = 0;
    for(uint32_t i = 0; i < statKeys.size(); ++i){
        ss.clear();
        ss.str("");
        ss << statMap[statKeys[i]];
        colStr = ss.str();
        colStr = util::replace(colStr, "\"", "");
        worksheet_write_string(sheet, row, 0, statKeys[i].c_str(), NULL);
        worksheet_write_number(sheet, row, 1, statMap[statKeys[i]], NULL);
        worksheet_write_number(sheet, row++, 2, (double)statMap[statKeys[i]]/statMap["TotalReads"], NULL);
        maxc1len = std::max(maxc1len, statKeys[i].length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
    worksheet_set_column(sheet, 2, 2, 6, NULL);
}

void genrpt::genBamQCSheet(lxw_worksheet* sheet, jsn::json& j){
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
        colStr = util::replace(colStr, "\"", "");
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

void genrpt::gen3UTRLenSheet(lxw_worksheet* sheet, jsn::json& j){
    size_t maxc1len = 0, maxc2len = 0;
    int row = 0;
    std::stringstream ss;
    std::string colStr;
    for(auto& iter: j.items()){
        ss.clear();
        ss.str("");
        if(iter.value().find("X1") == iter.value().end()){
            ss << 0;
        }else{
            ss << iter.value().at("X1");
        }
        colStr = ss.str();
        worksheet_write_string(sheet, row, 0, iter.key().c_str(), NULL);
        worksheet_write_number(sheet, row++, 1, std::atoi(colStr.c_str()), NULL);
        maxc1len = std::max(maxc1len, iter.key().length());
        maxc2len = std::max(maxc2len, colStr.length());
    }
    worksheet_set_column(sheet, 0, 0, maxc1len, NULL);
    worksheet_set_column(sheet, 1, 1, maxc2len, NULL);
}

void genrpt::genFinalReport(const genrpt::Stats& s, const genrpt::Options& opt){
    lxw_workbook* workbook = new_workbook(opt.outReportFile.c_str());
    std::ifstream fr;
    // bamqc sheets
    jsn::json jbamqc;
    fr.open(opt.bamqcResult.c_str());
    fr >> jbamqc;
    // sheet DupIncludeQC
    lxw_worksheet* sheet = workbook_add_worksheet(workbook, "DupIncludeQC");
    genBamQCSheet(sheet, jbamqc["DupIncludedQC"]);
    // sheet DupExcludeQC
    sheet = workbook_add_worksheet(workbook, "DupExcludeQC");
    genBamQCSheet(sheet, jbamqc["DupExcludeQC"]);
    fr.close();
    // sheet AllGeneExp
    sheet = workbook_add_worksheet(workbook, "AllGeneExp");
    int row = 0, col = 0;
    size_t maxc1len = 0, maxc2len = 0;
    std::stringstream oss;
    std::string prefix, suffix;
    worksheet_write_string(sheet, row, 0, "Gene", NULL);
    worksheet_write_string(sheet, row++, 1, "TPM", NULL);
    for(auto& e: s.geneExpress){
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
    fr.close();
    // sheet RnaClass
    fr.open(opt.filterResult.c_str());
    jsn::json jfilter;
    fr >> jfilter;
    sheet = workbook_add_worksheet(workbook, "RnaClass");
    genRNAClassSheet(sheet, jfilter["FilterResult"], jbamqc["DupIncludedQC"]);
    fr.close();
    // sheet 3UTRLen
    sheet = workbook_add_worksheet(workbook, "3UTRLen");
    gen3UTRLenSheet(sheet, jbamqc["IDP3UTRLen"]);
    // sheet Extra
    row = 0, col = 0;
    sheet = workbook_add_worksheet(workbook, "Extra");
    // DNA Pollution
    worksheet_write_string(sheet, row, col++, "DNAPolluteRate", NULL);
    double mappedReads = jbamqc["DupIncludedQC"]["EffectiveReads"];
    double filterGotReads = jbamqc["DupIncludedQC"]["TotalReads"];
    int32_t totReads = jfilter["FilterResult"]["Summary"]["ReadsIn"];
    worksheet_write_number(sheet, row++, col, 0.5 * (filterGotReads - mappedReads)/totReads, NULL);
    // Dup Rate
    col = 0;
    worksheet_write_string(sheet, row, col++, "DupRate", NULL);
    int32_t idpTTReads = jbamqc["DupIncludedQC"]["TotalReads"]; 
    int32_t ddpTTReads = jbamqc["DupExcludeQC"]["TotalReads"];
    worksheet_write_number(sheet, row++, col, (double)(idpTTReads - ddpTTReads)/idpTTReads, NULL);
    worksheet_set_column(sheet, 0, 0, 16, NULL);
    worksheet_set_column(sheet, 1, 1, 6, NULL);
    workbook_close(workbook);
    // plot graph
    std::string pltCMD = "integrity.py " + opt.bamqcResult + " " + util::joinpath(util::dirname(opt.outReportFile), util::replace(util::basename(opt.outReportFile), ".xlsx", ".png"));
    std::system(pltCMD.c_str());
}

#include "htslib/faidx.h"
#include "filereader.h"
#include "filewriter.h"
#include <iostream>
#include <cassert>
#include <sstream>
#include <map>

double gcContent(char* seq, int len){
    int gc = 0;
    for(int i = 0; i < len; ++i){
        if(seq[i] == 'G' || seq[i] == 'C' || seq[i] == 'g' || seq[i] == 'c'){
            ++gc;
        }
    }
    return double(gc)/double(len);
}

double avgDepth(std::vector<int>& depv, int beg, int end){
    int totDept = 0;
    for(int i = beg; i < end; ++i){
        totDept += depv[i];
    }
    return double(totDept)/double(end - beg);
}

void depMap(const char* file, std::map<std::string, std::vector<int>>& dmap){
    util::FileReader dr(file);
    std::string tstr;
    std::vector<std::string> vstr;
    while(dr.getline(tstr)){
        util::split(tstr, vstr, "\t");
        if(dmap.find(vstr[0]) == dmap.end()){
            dmap[vstr[0]] = {std::atoi(vstr[1].c_str())};
        }else{
            dmap[vstr[0]].push_back(std::atoi(vstr[1].c_str()));
        }
    }
}

int main(int argc, char** argv){
    if(argc < 6){
        std::cout << argv[0] << " <infa> <inbed> <step> <ddp.depth> <idp.depth> " << std::endl;
        return 0;
    }
    std::map<std::string, std::vector<int>> ddp;
    depMap(argv[4], ddp);
    std::map<std::string, std::vector<int>> idp;
    depMap(argv[5], idp);
    char* infa = argv[1];
    char* inreg = argv[2];
    int step = std::atoi(argv[3]);

    faidx_t* fai = fai_load(infa);
    util::FileReader fr(inreg);
    std::string line;
    std::vector<std::string> vs;
    int len;
    std::stringstream result;
    while(fr.getline(line)){
        util::split(line, vs, "\t");
        std::string name = vs[0];
        int beg = std::atoi(vs[1].c_str());
        int end = std::atoi(vs[2].c_str());
        char* s = faidx_fetch_seq(fai, name.c_str(), beg, end, &len);
        int count = 0;
        std::vector<double> gcv;
        std::vector<int> posv;
        std::vector<double> aidp;
        std::vector<double> addp;
        while(count + step <= len){
            gcv.push_back(gcContent(s + count, step));
            addp.push_back(avgDepth(ddp[name], count, count + step));
            aidp.push_back(avgDepth(idp[name], count, count + step));
            count += step;
            posv.push_back(count);
        }
        if(count < len){
            gcv.push_back(gcContent(s + count, len - count));
            addp.push_back(avgDepth(ddp[name], count, ddp[name].size()));
            aidp.push_back(avgDepth(idp[name], count, idp[name].size()));
            posv.push_back(len);
        }
        result << name;
        for(size_t i = 0; i < posv.size(); ++i){
            result << "\t" << posv[i];
        }
        result << "\nGCRatio";
        for(size_t i = 0; i < gcv.size(); ++i){
            result << "\t" << gcv[i];
        }
        result << "\nDDPMean";
        for(size_t i = 0; i < addp.size(); ++i){
            result << "\t" << addp[i];
        }
        result << "\nIDPMean";
        for(size_t i = 0; i < aidp.size(); ++i){
            result << "\t" << aidp[i];
        }
        result << "\n";
    }
    std::cout << result.rdbuf();
}

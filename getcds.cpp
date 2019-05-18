#include <set>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>
#include "util.h"
#include "filereader.h"
#include "filewriter.h"

struct ExonRecord{
    std::vector<std::pair<int, int>> exoncoord;
    int cdsStart = 0;
    int cdsEnd = 0;
    int trsStart = 0;
    int trsEnd = 0;
    int cdsLen = 0;
    int geneLen = 0;
    char strand = '+';
    int utr3Len = 0;
    int utr5Len = 0;
    int firstCDSIndex = 0;
    int lastCDSIndex = 0;
    std::string chromeName;
    std::string geneName;
    std::string trsName;
    std::string trsVer;

    ExonRecord() = default;
    ~ExonRecord() = default;
    void truncateByCDS(){
        for(uint32_t i = 0; i < exoncoord.size(); ++i){
            geneLen += exoncoord[i].second - exoncoord[i].first;
            if(cdsStart > exoncoord[i].second){
                utr5Len += exoncoord[i].second - exoncoord[i].first;
            }
            if(exoncoord[i].first <= cdsStart && exoncoord[i].second >= cdsStart){
                utr5Len += cdsStart - exoncoord[i].first;
                exoncoord[i].first = cdsStart;
                firstCDSIndex = i;
            }
            if(exoncoord[i].first <= cdsEnd && exoncoord[i].second >= cdsEnd){
                utr3Len += exoncoord[i].second - cdsEnd;
                exoncoord[i].second = cdsEnd;
                lastCDSIndex = i;
                break;
            }
        }
        for(uint32_t j = lastCDSIndex + 1; j < exoncoord.size(); ++j){
            geneLen += exoncoord[j].second - exoncoord[j].first;
            utr3Len += exoncoord[j].second - exoncoord[j].first;
        }
        cdsLen = geneLen - utr3Len - utr5Len;
        if(strand == '-'){
            std::swap(utr3Len, utr5Len);
        }
    }
    friend std::ostream& operator<<(std::ostream& os, const ExonRecord& rec){
        if(rec.exoncoord.size() == 0){
            return os;
        }
        os << "#" << rec.trsName << "\t" << rec.geneLen << "\t" << rec.cdsLen << "\t" << rec.utr5Len << "\t" << rec.utr3Len << "\n";
        os << ">" << rec.geneName << "_" << rec.trsName << "." << rec.trsVer << ",";
        os << rec.chromeName << ":" << rec.cdsStart << "-" << rec.cdsEnd << "\n";
        int i = 0;
        if(rec.strand == '+'){
            for(int k = rec.firstCDSIndex; k <= rec.lastCDSIndex; ++k){
                os << (++i) << "," << rec.exoncoord[k].first << "," << rec.exoncoord[k].second << "\n";
            }
            os << std::endl;
        }else{
            for(int k = rec.lastCDSIndex; k >= rec.firstCDSIndex; --k){
                os << (++i) << "," << rec.exoncoord[k].first << "," << rec.exoncoord[k].second << "\n";
            }
        }
        return os;
    }
};

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << std::string(argv[0]) << " <refgenewithaccs> <transcriptlist> <output> " << std::endl;
        std::exit(0);
    }
    util::FileReader freader = {std::string(argv[1])};
    util::FileWriter fwriter = {std::string(argv[3])};
    std::vector<std::string> vtrs;
    util::makeListFromFileByLine(std::string(argv[2]), vtrs);
    std::set<std::string> strs(vtrs.begin(), vtrs.end());
    std::string line;
    freader.getline(line);
    std::vector<std::string> vec;
    std::ostringstream oss;
    while(freader.getline(line)){
        util::split(line, vec, "\t");
        if(strs.find(vec[1]) == strs.end()){
            continue;
        }
        ExonRecord rec;
        rec.trsName = vec[1];
        rec.trsStart = std::atoi(vec[4].c_str());
        rec.trsEnd = std::atoi(vec[5].c_str());
        rec.cdsStart = std::atoi(vec[6].c_str());
        rec.cdsEnd = std::atoi(vec[7].c_str());
        rec.strand = vec[3][0];
        rec.chromeName = vec[2];
        rec.geneName = vec[12];
        rec.trsVer = vec[16];
        std::vector<std::string> exonStart;
        std::vector<std::string> exonEnd;
        util::split(vec[9], exonStart, ",");
        util::split(vec[10], exonEnd, ",");
        std::vector<int> startCoord;
        std::vector<int> endCoord;
        util::strvec2intvec(exonStart, startCoord);
        util::strvec2intvec(exonEnd, endCoord);
        util::vec2pairvec(startCoord, endCoord, rec.exoncoord);
        rec.truncateByCDS();
        oss.clear();
        oss.str("");
        oss << rec;
        fwriter.writeString(oss.str());
    }
}


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
    std::vector<std::pair<size_t, size_t>> exoncoord;
    size_t cdsStart;
    size_t cdsEnd;
    int firstCDSIndex;
    int lastCDSIndex;
    std::string chromeName;
    std::string geneName;
    std::string trsName;
    std::string trsVer;

    ExonRecord() = default;
    ~ExonRecord() = default;
    void truncateByCDS(){
        firstCDSIndex = 0;
        lastCDSIndex = 0;
        for(int i = 0; i < exoncoord.size(); ++i){
            if(exoncoord[i].first <= cdsStart && exoncoord[i].second >= cdsStart){
                exoncoord[i].first = cdsStart;
                firstCDSIndex = i;
            }
            if(exoncoord[i].first <= cdsEnd && exoncoord[i].second >= cdsEnd){
                exoncoord[i].second = cdsEnd;
                lastCDSIndex = i;
                break;
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const ExonRecord& rec){
        if(rec.exoncoord.size() == 0){
            return os;
        }
        os << ">" << rec.geneName << "_" << rec.trsName << "." << rec.trsVer << ",";
        os << rec.chromeName << ":" << rec.cdsStart << "-" << rec.cdsEnd << "\n";
        int i = 0;
        for(int k = rec.firstCDSIndex; k <= rec.lastCDSIndex; ++k){
            os << (++i) << "," << rec.exoncoord[k].first << "," << rec.exoncoord[k].second << "\n";
        }
        os << std::endl;
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
        rec.cdsStart = std::atoi(vec[6].c_str());
        rec.cdsEnd = std::atoi(vec[7].c_str());
        rec.chromeName = vec[2];
        rec.geneName = vec[12];
        rec.trsVer = vec[16];
        std::vector<std::string> exonStart;
        std::vector<std::string> exonEnd;
        util::split(vec[9], exonStart, ",");
        util::split(vec[10], exonEnd, ",");
        std::vector<size_t> startCoord;
        std::vector<size_t> endCoord;
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


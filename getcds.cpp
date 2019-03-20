#include <set>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>
#include "util.h"
#include "filereader.h"
#include "filewriter.h"

/** object to store some infomation ofan record of refgene with transcript number added */
struct ExonRecord{
    std::vector<std::pair<size_t, size_t>> exoncoord; ///< exon coordinates pair records in genome
    size_t cdsStart;          ///< CDS starting coordinates in genome
    size_t cdsEnd;            ///< CDS ending coordinates in genome
    int firstCDSIndex;        ///< first CDS index in exoncoord
    int lastCDSIndex;         ///< last CDS index in exoncorod
    std::string chromeName;   ///< chromesome name of this gene                
    std::string geneName;     ///< gene name 
    std::string trsName;      ///< transcript name
    std::string trsVer;       ///< transcript version

    /** construct an ExonRecord object */
    ExonRecord() = default;
    /** destroy an ExonRecord object */
    ~ExonRecord() = default;
   
    /** get the cds index in exoncoord and adjust the edging coordinates */
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

    /** output an exon record to ostream
     * @param os reference to a ostream
     * @param rec reference to a ExonRecord object
     * @return reference to a ostream
     */
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

/** main function to accept arguments
 * @argc number of arguments provided externally
 * @argv array of char* store the externally provided arguments
 * @return 0 if exit successfuly
 */
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

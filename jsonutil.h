#ifndef JSON_UTIL_H
#define JSON_UTIL_H

#include <fstream>
#include <string>

/** json output utilities */
namespace jsonutil{
    /** output a key-val record to ofstream in json format
     * @param ofs std::ofstream to output to
     * @param padding padding before each record
     * @param key the key of a json record
     * @param val the value of a json record
     */
    template<typename T>
    inline void writeRecord(std::ofstream& ofs, const std::string& padding, const std::string& key, const T& val){
        ofs << padding << "\t\"" << key << "\": " << val << ",\n";
    }

    /** get current system time
     * @return current system time
     */
    inline std::string getCurrentSystemTime(){
        time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        struct tm* ptm = std::localtime(&tt);
        char date[60] = {0};
        std::sprintf(date, "%d-%02d-%02d %02d:%02d:%02d",
                (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
                (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
        return std::string(date);
    }
}

#endif

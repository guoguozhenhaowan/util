#ifndef HTML_UTIL_H
#define HTML_UTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

/** some useful utilities for html output */
namespace htmlutil{
    /** output a row with two columns as a html table row to ofstream
     * @param ofs std::ofstream to output to
     * @param key first field of a row
     * @param val second field of a row
     */
    template<typename T>
    void outputTableRow(std::ofstream& ofs, const std::string& key, T val){
        std::stringstream ss;
        ss << val;
        ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + ss.str() + "</td></tr>\n";
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

    /** format a big number with KMGTP units
     * @param number a non-negative number to be formatted
     * @return formatted number string
     */
    inline std::string formatNumber(const size_t& number){
        double num = (double)number;
        const std::string unit[6] = {"", "K", "M", "G", "T", "P"};
        int order = 0;
        while(number > 1000.0){
            order += 1;
            num /= 1000.0;
        }
        return std::to_string(num) + " " + unit[order];
    }

    /** get percents number of nuerator/denominator
     * @param numerator numerator
     * @param denominator denominator
     * @return percents number string
     */
    template<typename T>
    std::string getPercentsStr(T numerator, T denominator){
        if(denominator == 0){
            return "0.0";
        }
        return std::to_string((double)numerator * 100.0 / (double) denominator);
    }
}

#endif

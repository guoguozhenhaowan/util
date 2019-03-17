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
}

#endif

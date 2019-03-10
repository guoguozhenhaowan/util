#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

namespace util{

    /** whether a string starts with some substring
     * @param str whole string
     * @param pre substring
     */
    inline bool starts_with(const std::string& str, const std::string& pre){
        if(str.length() < pre.length()){
            return false;
        }else{
            return std::equal(pre.begin(), pre.end(), str.begin());
        }
    }

    /** whether a string ends with some substring
     * @param str whole string
     * @param suf substring
     */
    inline bool ends_with(const std::string& str, const std::string& suf){
        if(str.length() < suf.length()){
            return false;
        }else{
            return std::equal(suf.rbegin(), suf.rend(), str.rbegin());
        }
    }

    /** get rid of the leading and ending white space characters of a string
     * @param str string to be stripped in both ends
     */
    inline std::string strip(const std::string& str){
        std::string::size_type ipos = str.find_first_not_of(" \t\n\v\f\r");
        if(ipos == std::string::npos){
            return "";
        }
        std::string::size_type epos = str.find_last_not_of(" \t\n\v\f\r");
        if(epos == ipos){
            return str.substr(ipos);
        }else{
            return str.substr(ipos, epos - ipos + 1);
        }
    }

    /** get rid of the leading white space characters of a string
     * @param str string to be stripped from front
     */
    inline std::string lstrip(const std::string& str){
        std::string::size_type pos = str.find_first_not_of(" \t\n\v\f\r");
        if(pos == std::string::npos){
            return "";
        }
        return str.substr(pos);
    }

    /** get rid of the trailling white space characters of a string
     * @param str string to be stripped from back
     */
    inline std::string rstrip(const std::string& str){
        std::string::size_type pos = str.find_last_not_of(" \t\n\v\f\r");
        if(pos == std::string::npos){
           return "";
        }
        return str.substr(0, pos + 1);
    } 

    /** split a string by predefined seperator into a vector
     * @param str string
     * @param vec vector to store the split results
     * @param sep seperators, can contain a series of seperators
     */
    inline void split(const std::string& str, std::vector<std::string>& vec, std::string sep = " "){
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

    /** replace a substr apearing in a string with another string
     * @param str string 
     * @param pat substr of string to be replaced
     * @param des string to be used to replaced with pat
     */
    inline std::string replace(const std::string&str, const std::string& pat, const std::string& des){
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

#endif

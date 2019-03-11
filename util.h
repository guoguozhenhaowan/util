#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <cerrno>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sys/stat.h>

namespace util{

    /** whether a string starts with some substring
     * @param str whole string
     * @param pre substring
     * @return true if str starts with pre
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
     * @return true if str ends with suf
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
     * @return a string with white spaces stripped in both ends
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

    /** get rid of the left leading white space characters of a string
     * @param str string to be stripped from front
     * @return a string with left leading white spaces stripped
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
     * @return a string with right ending white spaces stripped
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
     * @return a string with each pat replaced by des
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

    /** get the basename of a path string
     * @param path name
     * @return basename of path
     */
    inline std::string basename(const std::string& path){
        if(path.find_first_of(" \t\n\v\f\r") != std::string::npos){
            return "";
        }
        std::string::size_type pos1 = path.find_last_of("/\\");
        if(pos1 == std::string::npos){
            return path;
        }
        std::string::size_type pos2 = path.find_last_not_of("/\\");
        if(pos2 == path.size() - 1){
            return path.substr(pos1 + 1, pos2 - pos1);
        }
        std::string::size_type pos3 = path.find_last_of("/\\", pos2);
        if(pos3 == std::string::npos){
            return path.substr(0, pos2 + 1);
        }
        return path.substr(pos3 + 1, pos2 - pos3);
    }

    /** get the dirname of a path string
     * @param path name
     * @return dirname of path
     */
    inline std::string dirname(const std::string& path){
        std::string::size_type pos = path.find_last_of("/\\");
        if(pos == std::string::npos){
#ifdef _WIN32
            return ".\\";
#else
            return "./";
#endif
        }
        if(pos == path.size() - 1){
            std::string::size_type pos1 = path.find_last_not_of("/\\");
            if(pos1 == std::string::npos){
#ifdef _WIN32
                return "\\";
#else
                return "/";
#endif
            }
            std::string::size_type pos2 = path.find_last_of("/\\", pos1);
            if(pos2 == std::string::npos){
#ifdef _WIN32
                return ".\\";
#else
                return "./";
#endif
            }else{
                return path.substr(0, pos2 + 1);
            }
        }else{
            return path.substr(0, pos + 1);
        }
    }

    /** join dirname and basename into a path
     * @param dirname dirname of the path
     * @param basename basename of the path
     * @return full path
     */
    inline std::string joinpath(const std::string& dirname, const std::string& basename){
#ifdef _WIN32
        return dirname + "\\" + basename;
#else
        return dirname + "/" + basename;
#endif
    }

    /** check a string is a regular file or not
     * @param path string of a file/directory
     * @return true if path is an existing regular file
     */
    inline bool isfile(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFREG);
#else
        struct stat info;
        if(stat(path.c_str(), &info) !=0){
            return false;
        }
        return (info.st_mode & S_IFREG);
#endif
    }

    /** check a string is a directory or not
     * @param path string of a file/directory
     * @return true if path is an existing path
     */
    inline bool isdir(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFDIR);
#else
        struct stat info;
        if(stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & S_IFDIR);
#endif
    }
    
    /** check a file/directory exists or not
     * @param path string of a file/directory
     * @return true if path exists
     */
    inline bool exists(const std::string& path){
        return util::isdir(path) || util::isfile(path);
    }

    /** make directories recursively
     * @param path path of directory to be created
     * @return true if make directories successfully
     */
    bool makedir(const std::string& path){
#ifdef _WIN32
        int ret = _mkdir(path.c_str());
#else
        mode_t mode = 0755;
        int ret = mkdir(path.c_str(), mode);
#endif
        if(ret == 0){
            return true;
        }
        switch(errno){
            case ENOENT:
                {
                    std::cout << util::dirname(path) << std::endl;
                    if(!util::makedir(util::dirname(path))){
                        return false;
                    }
                }
#ifdef _WIN32
                return 0 == _mkdir(path.c_str());
#else
                return 0 == mkdir(path.c_str(), mode);
#endif
            case EEXIST:
                return util::isdir(path);
            default:
                return false;
        }
    }
}

#endif

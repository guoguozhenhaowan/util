#ifndef FILEWRITER_HPP
#define FILEWRITER_HPP

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <zlib.h>

/** Class to write to gz file or ofstream */
class FileWriter{
    std::string mFilename;  ///< output filename
    gzFile mGzFile;         ///< gzFile file handler
    std::ofstream* mStream; ///< pointer to ofstream
    bool mZipped;           ///< output file is mZipped or not
    int mCompressLevel;     ///< compression level for gz file
    bool mNeedClose;        ///< needed to be closed or not
    
    public:
    /** FileWriter constructor
     * @param filename output filename
     * @param compression compression level for gzFile
     */
    FileWriter(const std::string& filename, const int& compression = 3){
        mCompressLevel = compression;
        mFilename = filename;
        mGzFile = NULL;
        mZipped = false;
        mNeedClose = true;
        init();
    }
        
    /** FileWriter constructor
     * @param stream pointer to ofstream
     */
    FileWriter(std::ofstream* stream){
        mGzFile = NULL;
        mZipped = false;
        mStream = stream;
        mNeedClose = false; 
    }
    
    /** FileWriter constructor
     * @param gzfile gzFile handler
     */
    FileWriter(gzFile gzfile){
        mStream = NULL;
        mGzFile = gzfile;
        mZipped = true;
        mNeedClose = false;
    }
    
    /** FileWriter destructor */
    ~FileWriter(){
        if(mNeedClose){
            close();
        }
    }

    /** Whether the output file is mZipped or not
     * @return true if output filename ends with .gz
     */
    inline bool isZipped(){
        return mZipped;
    }
    
    /** write a string to file without append \\n
     * @param str string to be written 
     * @return true if successfully written
     */
    inline bool writeString(const std::string& str){
        const char* cstr = str.c_str();
        size_t size = str.length();
        size_t written = 0;
        bool status = true;
        if(mZipped){
            written = gzwrite(mGzFile, cstr, size);
            status = size == written;
        }else{
            mStream->write(cstr, size);
            status = !mStream->fail();
        }
        return status;
    }
    
    /** write a string to file with additional \\n appended 
     * @param linestr string to be written
     * @return true if successfully written
     */
    inline bool writeLine(const std::string& linestr){
        const char* line = linestr.c_str();
        size_t size = linestr.length();
        size_t written = 0;
        bool status = true;
        if(mZipped){
            written = gzwrite(mGzFile, line, size);
            gzputc(mGzFile, '\n');
            status = size == written;
        }else{
            mStream->write(line, size);
            mStream->put('\n');
            status = !mStream->fail();
        }
        return status;
    }
    
    /** write a C string to file
     * @param cstr pointer to char(C string)
     * @param size length of the C string to be written
     * @return true if successfully written
     */
    inline bool write(char* cstr, size_t size){
        size_t written = 0;
        bool status = true;
        if(mZipped){
            written = gzwrite(mGzFile, cstr, size);
            status = size == written;
        }else{
            mStream->write(cstr, size);
            status = !mStream->fail();
        }
        return status;
    }
    
    /** get filename of output file
     * @return output filename
     */
    inline std::string getFilename(){
        return mFilename;
    }

    /** initialize FileWriter, detect file format and open file handler
     */
    inline void init(){
        if(util::ends_with(mFilename, ".gz")){
            mGzFile = gzopen(mFilename.c_str(), "w");
            gzsetparams(mGzFile, mCompressLevel, Z_DEFAULT_STRATEGY);
            gzbuffer(mGzFile, 1024 * 1024);
            mZipped = true;
        }else{
            mStream = new std::ofstream();
            mStream->open(mFilename.c_str(), std::ios::out);
            mZipped = false;
        }
    }
    
    /** flush buffer and close file handler
     */
    inline void close(){
        if(mZipped){
            if(mGzFile){
                gzflush(mGzFile, Z_FINISH);
                gzclose(mGzFile);
                mGzFile = NULL;
            }
        }else if(mStream){
            if(mStream->is_open()){
                mStream->flush();
                mStream = NULL;
            }
        }
    }
};

#endif

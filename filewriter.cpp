#include "filewriter.h"

namespace util{
    FileWriter::FileWriter(const std::string& filename, const int& compression){
        mCompressLevel = compression;
        mFilename = filename;
        mGzFile = NULL;
        mZipped = false;
        mNeedClose = true;
        init();
    }
    
    FileWriter::FileWriter(std::ofstream* stream){
        mGzFile = NULL;
        mZipped = false;
        mStream = stream;
        mNeedClose = false;
    }
    
    FileWriter::FileWriter(gzFile gzfile){
        mStream = NULL;
        mGzFile = gzfile;
        mZipped = true;
        mNeedClose = false;
    }
    
    FileWriter::~FileWriter(){
        if(mNeedClose){
            close();
        }
    }
    
    std::string FileWriter::getFilename(){
        return mFilename;
    }
    
    void FileWriter::init(){
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
    
    bool FileWriter::writeLine(const std::string& linestr){
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
    
    bool FileWriter::writeString(const std::string& str){
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
    
    bool FileWriter::write(char* cstr, size_t size){
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
    
    void FileWriter::close(){
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
    
    bool FileWriter::isZipped(){
        return mZipped;
    }
}

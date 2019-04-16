#include "util.h"
#include "filereader.h"
#include <cstring>

namespace util{
    FileReader::FileReader(const std::string& filename){
        mFileName = filename;
        mGzipFile = NULL;
        mFile = NULL;
        mStdinMode = false;
        mReadBufSize = (1 << 20);
        mBuf = new char[mReadBufSize];
        mBufDataLen = 0;
        mBufUsedLen = 0;
        mNoLineBreakAtEnd = false;
        init();
    }
    
    FileReader::~FileReader(){
        close();
        delete mBuf;
        mBuf = nullptr;
    }
    
    bool FileReader::hasNoLineBreakAtEnd(){
        return mNoLineBreakAtEnd;
    }
    
    void FileReader::readToBuf(){
        if(isZipped()){
            mBufDataLen = gzread(mGzipFile, mBuf, mReadBufSize);
            if(mBufDataLen == -1){
                std::cerr << "Error to read gzip file" << std::endl;
            }
        }else{
            mBufDataLen = std::fread(mBuf, 1, mReadBufSize, mFile);
        }
        mBufUsedLen = 0;
        if(mBufDataLen < mReadBufSize){
            if(mBuf[mBufDataLen - 1] != '\n'){
                mNoLineBreakAtEnd = true;
            }
        }
    }
    
    void FileReader::init(){
        if(util::ends_with(mFileName, ".gz")){
            mGzipFile = ::gzopen(mFileName.c_str(), "r");
            mZipped = true;
            ::gzrewind(mGzipFile);
        }else{
            if(mFileName == "/dev/stdin"){
                mFile = stdin;
            }else{
                mFile = std::fopen(mFileName.c_str(), "rb");
            }
            if(mFile == NULL){
                util::error_exit("Failed to open file: " + mFileName);
            }
            mZipped = false;
        }
        readToBuf();
    }
    
    void FileReader::getBytes(size_t& bytesRead, size_t& bytesTotal){
        if(mZipped){
            bytesRead = ::gzoffset(mGzipFile);
        }else{
            bytesRead = std::ftell(mFile);
        }
        
        // use another ifstream without affecting the current reader
        std::ifstream is(mFileName);
        is.seekg(0, is.end);
        bytesTotal = is.tellg();
    }
    
    void FileReader::clearLineBreaks(char* line){
        // trim \n, \r or \r\n in the tail
        int len = std::strlen(line);
        if(len >= 2){
            if(line[len - 1] == '\n' || line[len - 1] == '\r'){
                line[len - 1] = '\0';
                if(line[len - 2] == '\r'){
                    line[len - 2] = '\0';
                }
            }
        }
    }

    bool FileReader::getline(std::string& line){
        if(mBufUsedLen >= mBufDataLen && eof()){
            return false;
        }
        if(mZipped && mGzipFile == NULL){
            return false;
        }
        line = getlineFromBuffer();
        return true;
    }

    std::string FileReader::getlineFromBuffer(){
        static int c = 0;
        ++c;
        int copied = 0;
        int start = mBufUsedLen;
        int end = start;
        
        // look for '\r' or '\n' until the end of mBuf
        while(end < mBufDataLen){
            if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                ++end;
            }else{
                break;
            }
        }
    
        // if '\r' or '\n' found(this line well contained in this mBuf)
        // or this is the last mBuf of file
        if(end < mBufDataLen || mBufDataLen < mReadBufSize){
            int len = end - start;
            std::string line(mBuf + start, len);
            ++end;
            if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                ++end;
            }
            mBufUsedLen = end;
            return line;
        }
    
        // if '\r' or '\n' not found && this is not the last mBuf of file
        // then this line is not contained in this mBuf, we should read new mBuf
        std::string str(mBuf + start, mBufDataLen - start);
        while(true){
            readToBuf();
            start = 0;
            end = 0;
            // look for '\r' or '\n' until the end of mBuf
            while(end < mBufDataLen){
                if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                    ++end;
                }else{
                    break;
                }
            }
    
            // if '\r' or '\n' found(this line well contained in this mBuf)
            // or this is the last mBuf of file
            if(end < mBufDataLen || mBufDataLen < mReadBufSize){
                int len = end - start;
                str.append(mBuf + start, len);
                ++end;
                if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                    ++end;
                }
                mBufUsedLen = end;
                return str;
            }
    
            // if '\r' or '\n' not found && this is not the last mBuf of file
            // then this line is not contained in this mBuf, we should read new mBuf
            str.append(mBuf + start, mBufDataLen);
        }
        return std::string();
    }
    
    bool FileReader::eof(){
        if(mZipped){
            return ::gzeof(mGzipFile);
        }else{
            return std::feof(mFile);
        }
    }
    
    void FileReader::close(){
        if(mZipped && mGzipFile){
            ::gzclose(mGzipFile);
            mGzipFile = NULL;
        }else if(mFile){
            std::fclose(mFile);
            mFile = NULL;
        }else{
            return;
        }
    }
    
    
    bool FileReader::isZippedFile(const std::string& filename){
        // just use the suffix of filename to determine its format
        if(util::ends_with(filename, ".fastq.gz") || util::ends_with(filename, ".fq.gz") ||
           util::ends_with(filename, ".fasta.gz") || util::ends_with(filename, ".fa.gz")){
            return true;
        }
        return false;
    }
    
    bool FileReader::isZipped(){
        return mZipped;
    }
}

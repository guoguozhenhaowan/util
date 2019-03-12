#include "util.h"
#include "fqreader.h"
#include <cstring>

#define FQ_BUF_SIZE (1<<20)

FqReader::FqReader(const std::string& fileName, const bool& hasQuality, const bool& phread64){
    this->fileName = fileName;
    this->gzipFile = NULL;
    this->pfile = NULL;
    this->stdinMode = false;
    this->phread64 = phread64;
    this->hasQuality = hasQuality;
    this->buf = new char[FQ_BUF_SIZE];
    this->bufDataLen = 0;
    this->bufUsedLen = 0;
    this->noLineBreakAtEnd = false;
    this->init();
}

FqReader::~FqReader(){
    this->close();
    delete this->buf;
    this->buf = nullptr;
}

bool FqReader::hasNoLineBreakAtEnd(){
    return this->noLineBreakAtEnd;
}

void FqReader::readToBuf(){
    if(this->isZipped()){
        this->bufDataLen = ::gzread(this->gzipFile, this->buf, FQ_BUF_SIZE);
        if(this->bufDataLen == -1){
            std::cerr << "Error to read gzip file" << std::endl;
        }
    }else{
        this->bufDataLen = std::fread(this->buf, 1, FQ_BUF_SIZE, this->pfile);
    }
    this->bufUsedLen = 0;
    if(this->bufDataLen < FQ_BUF_SIZE){
        if(this->buf[bufDataLen - 1] != '\n'){
            this->noLineBreakAtEnd = true;
        }
    }
}

void FqReader::init(){
    if(util::ends_with(this->fileName, ".gz")){
        this->gzipFile = ::gzopen(this->fileName.c_str(), "r");
        this->zipped = true;
        ::gzrewind(gzipFile);
    }else{
        if(this->fileName == "/dev/stdin"){
            this->pfile = stdin;
        }else{
            this->pfile = std::fopen(this->fileName.c_str(), "rb");
        }
        if(this->pfile == NULL){
            util::error_exit("Failed to open file: " + this->fileName);
        }
        this->zipped = false;
    }
    this->readToBuf();
}

void FqReader::getBytes(size_t& bytesRead, size_t& bytesTotal){
    if(this->zipped){
        bytesRead = ::gzoffset(this->gzipFile);
    }else{
        bytesRead = std::ftell(this->pfile);
    }
    
    // use another ifstream without affecting the current reader
    std::ifstream is(this->fileName);
    is.seekg(0, is.end);
    bytesTotal = is.tellg();
}

void FqReader::clearLineBreaks(char* line){
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

std::string FqReader::getLine(){
    static int c = 0;
    ++c;
    int copied = 0;
    int start = this->bufUsedLen;
    int end = start;
    
    // look for '\r' or '\n' until the end of buf
    while(end < this->bufDataLen){
        if(this->buf[end] != '\r' && this->buf[end] != '\n'){
            ++end;
        }else{
            break;
        }
    }

    // if '\r' or '\n' found(this line well contained in this buf)
    // or this is the last buf of file
    if(end < this->bufDataLen || this->bufDataLen < FQ_BUF_SIZE){
        int len = end - start;
        std::string line(this->buf + start, len);
        ++end;
        if(end < this->bufDataLen - 1 && this->buf[end - 1] == '\r' && this->buf[end] == '\n'){
            ++end;
        }
        this->bufUsedLen = end;
        return line;
    }

    // if '\r' or '\n' not found && this is not the last buf of file
    // then this line is not contained in this buf, we should read new buf
    std::string str(this->buf + start, bufDataLen - start);
    while(true){
        this->readToBuf();
        start = 0;
        end = 0;
        // look for '\r' or '\n' until the end of buf
        while(end < this->bufDataLen){
            if(this->buf[end] != '\r' && this->buf[end] != '\n'){
                ++end;
            }else{
                break;
            }
        }

        // if '\r' or '\n' found(this line well contained in this buf)
        // or this is the last buf of file
        if(end < this->bufDataLen || this->bufDataLen < FQ_BUF_SIZE){
            int len = end - start;
            str.append(this->buf + start, len);
            ++end;
            if(end < this->bufDataLen - 1 && this->buf[end - 1] == '\r' && this->buf[end] == '\n'){
                ++end;
            }
            this->bufUsedLen = end;
            return str;
        }

        // if '\r' or '\n' not found && this is not the last buf of file
        // then this line is not contained in this buf, we should read new buf
        str.append(this->buf + start, this->bufDataLen);
    }
    return std::string();
}

bool FqReader::eof(){
    if(this->zipped){
        return ::gzeof(this->gzipFile);
    }else{
        return std::feof(this->pfile);
    }
}

Read* FqReader::read(){
    if(this->zipped && this->gzipFile == NULL){
        return NULL;
    }
    if(this->bufUsedLen >= this->bufDataLen && this->eof()){
        return NULL;
    }

    std::string name = this->getLine();
    while((name.empty() && !(this->bufUsedLen >= this->bufDataLen && this->eof())) || (!name.empty() && name[0] !='@')){
        name = this->getLine();
    }
    if(name.empty()){
        return NULL;
    }

    std::string sequence = this->getLine();
    std::string strand = this->getLine();
    // some fq has no quality, then construct the quality string with all 'K'
    if(!this->hasQuality){
        std::string quality = std::string(sequence.length(), 'K');
        return new Read(name, sequence, strand, quality, this->phread64);
    }else{
        std::string quality = this->getLine();
        if(quality.length() != sequence.length()){
            std::cerr << "Error: base sequnce and quality sequence have different length: \n";
            std::cerr << name << "\n";
            std::cerr << sequence << "\n";
            std::cerr << quality << "\n";
            std::cerr << strand << "\n";
            return NULL;
        }
        return new Read(name, sequence, strand, quality, this->phread64);
    }
    return NULL;
}

void FqReader::close(){
    if(this->zipped && this->gzipFile){
        ::gzclose(this->gzipFile);
        this->gzipFile = NULL;
    }else if(this->pfile){
        std::fclose(this->pfile);
        this->pfile = NULL;
    }else{
        return;
    }
}


bool FqReader::isZippedFq(const std::string& filename){
    // just use the suffix of filename to determine its format
    if(util::ends_with(filename, ".fastq.gz") || util::ends_with(filename, ".fq.gz") ||
       util::ends_with(filename, ".fasta.gz") || util::ends_with(filename, ".fa.gz")){
        return true;
    }
    return false;
}

bool FqReader::isfq(const std::string& filename){
    // just use the suffix of filename to determine its format
    if(util::ends_with(filename, ".fastq") || util::ends_with(filename, ".fq") ||
       util::ends_with(filename, ".fasta") || util::ends_with(filename, ".fa")){
        return true;
    }
    return false;
}

bool FqReader::isZipped(){
    return this->zipped;
}

FqReaderPair::FqReaderPair(const std::string& lname, const std::string& rname, const bool& hasQual,
                           const bool& phread64, const bool& interleaved){
    this->interleaved = interleaved;
    this->left = new FqReader(lname, hasQual, phread64);
    if(interleaved){
        this->right = NULL;
    }else{
        this->right = new FqReader(rname, hasQual, phread64);
    }
}

FqReaderPair::~FqReaderPair(){
    if(this->left){
        delete this->left;
        this->left = NULL;
    }
    if(this->right){
        delete this->right;
        this->right = NULL;
    }
}

ReadPair* FqReaderPair::read(){
    Read* pr1 = this->left->read();
    Read* pr2 = NULL;
    if(this->interleaved){
        pr2 = this->left->read();
    }else{
        pr2 = this->right->read();
    }
    if(!pr1 || !pr2){
        return NULL;
    }else{
        return new ReadPair(pr1, pr2);
    }
}

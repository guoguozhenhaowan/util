#include "writer.h"

Writer::Writer(const std::string& filename, const int& compression){
    this->compressLevel = compression;
    this->filename = filename;
    this->zipfile = NULL;
    this->zipped = false;
    this->needClose = true;
    this->init();
}

Writer::Writer(std::ofstream* ofs){
    this->zipfile = NULL;
    this->zipped = false;
    this->ofs = ofs;
    this->needClose = false;
}

Writer::Writer(gzFile gzfile){
    this->ofs = NULL;
    this->zipfile = gzfile;
    this->zipped = true;
    this->needClose = false;
}

Writer::~Writer(){
    if(needClose){
        this->close();
    }
}

std::string Writer::getFilename(){
    return this->filename;
}

void Writer::init(){
    if(util::ends_with(this->filename, ".gz")){
        this->zipfile = gzopen(this->filename.c_str(), "w");
        gzsetparams(this->zipfile, this->compressLevel, Z_DEFAULT_STRATEGY);
        gzbuffer(this->zipfile, 1024 * 1024);
        this->zipped = true;
    }else{
        this->ofs = new std::ofstream();
        this->ofs->open(this->filename.c_str(), std::ios::out);
        this->zipped = false;
    }
}

bool Writer::writeLine(const std::string& linestr){
    const char* line = linestr.c_str();
    size_t size = linestr.length();
    size_t written = 0;
    bool status = true;
    if(this->zipped){
        written = gzwrite(this->zipfile, line, size);
        gzputc(this->zipfile, '\n');
        status = size == written;
    }else{
        this->ofs->write(line, size);
        this->ofs->put('\n');
        status = !this->ofs->fail();
    }
    return status;
}

bool Writer::writeString(const std::string& str){
    const char* cstr = str.c_str();
    size_t size = str.length();
    size_t written = 0;
    bool status = true;
    if(this->zipped){
        written = gzwrite(this->zipfile, cstr, size);
        status = size == written;
    }else{
        this->ofs->write(cstr, size);
        status = !this->ofs->fail();
    }
    return status;
}

bool Writer::write(char* cstr, size_t size){
    size_t written = 0;
    bool status = true;
    if(this->zipped){
        written = gzwrite(this->zipfile, cstr, size);
        status = size == written;
    }else{
        this->ofs->write(cstr, size);
        status = !this->ofs->fail();
    }
    return status;
}

void Writer::close(){
    if(this->zipped){
        if(this->zipfile){
            gzflush(this->zipfile, Z_FINISH);
            gzclose(this->zipfile);
            this->zipfile = NULL;
        }
    }else if(this->ofs){
        if(this->ofs->is_open()){
            this->ofs->flush();
            this->ofs = NULL;
        }
    }
}

bool Writer::isZipped(){
    return this->zipped;
}


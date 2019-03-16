#include "writerthread.h"

WriterThread::WriterThread(const std::string& filename, int compression){
    this->writer = NULL;
    this->inputCounter = 0;
    this->outputCounter = 0;
    this->inputCompleted = false;
    this->filename = filename;

    this->ringBuffer = new char*[compar::PACK_NUM_LIMIT];
    std::memset(this->ringBuffer, 0, sizeof(char*) * compar::PACK_NUM_LIMIT);
    this->ringBufferSizes = new size_t[compar::PACK_NUM_LIMIT];
    std::memset(this->ringBufferSizes, 0, sizeof(size_t) * compar::PACK_NUM_LIMIT);
    this->iniWriter(filename, compression);
}

WriterThread::~WriterThread(){
    this->cleanup();
    delete this->ringBuffer;
}

bool WriterThread::isCompleted(){
    return this->inputCompleted && (this->outputCounter == this->inputCounter);
}

bool WriterThread::setInputCompleted(){
    this->inputCompleted = true;
    return true;
}

void WriterThread::output(){
    if(this->outputCounter >= this->inputCounter){
        ::usleep(100);
    }
    while(this->outputCounter < this->inputCounter){
        this->writer->write(this->ringBuffer[this->outputCounter], 
                            this->ringBufferSizes[this->outputCounter]);
        delete this->ringBuffer[this->outputCounter];
        this->ringBuffer[this->outputCounter] = NULL;
        ++this->outputCounter;
    }
}

void WriterThread::input(char* cstr, size_t size){
    this->ringBuffer[this->inputCounter] = cstr;
    this->ringBufferSizes[this->inputCounter] = size;
    ++this->inputCounter;
}

void WriterThread::cleanup(){
    this->deleteWriter();
}

void WriterThread::deleteWriter(){
    if(this->writer != NULL){
        delete this->writer;
        this->writer = NULL;
    }
}

void WriterThread::iniWriter(const std::string& filename, int compression){
    this->deleteWriter();
    this->writer = new Writer(filename, compression);
}

void WriterThread::iniWriter(std::ofstream* ofs){
    this->deleteWriter();
    this->writer = new Writer(ofs);
}

void WriterThread::iniWriter(gzFile gzfile){
    this->deleteWriter();
    this->writer = new Writer(gzfile);
}

size_t WriterThread::bufferLength(){
    return this->inputCounter - this->outputCounter;
}

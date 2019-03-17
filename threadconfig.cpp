#include "threadconfig.h"

ThreadConfig::ThreadConfig(ThreadOpt* opt, FilterOpt* fopt, int threadId, bool paired){
    this->opt = opt;
    this->threadId = threadId;
    this->workingSplit = threadId;
    this->currentSplitReads = 0;
    this->preStats1 = new Stats(this->opt->estimatedReadLen);
    this->postStats1 = new Stats(this->opt->estimatedReadLen);
    if(paired){
        this->preStats2 = new Stats(this->opt->estimatedReadLen);
        this->postStats2 = new Stats(this->opt->estimatedReadLen);
    }else{
        this->preStats2 = NULL;
        this->postStats2 = NULL;
    }
    this->writer1 = NULL;
    this->writer2 = NULL;
    this->filterResult = new FilterResult(this->fopt, paired);
    this->canBeStopped = false;
}

ThreadConfig::~ThreadConfig(){
    this->cleanup();
}

void ThreadConfig::cleanup(){
    if(this->opt->enableSplit && this->opt->splitByFileNumber){
        this->writeEmptyFilesForSplitting();
    }
    this->deleteWriter();
}

void ThreadConfig::deleteWriter(){
    if(this->writer1 != NULL){
        delete this->writer1;
        this->writer1 = NULL;
    }
    if(this->writer2 != NULL){
        delete this->writer2;
        this->writer2 = NULL;
    }
}

void ThreadConfig::initWriter(const std::string& filename){
    this->deleteWriter();
    this->writer1 = new Writer(filename, this->opt->compression);
}

void ThreadConfig::initWriter(const std::string& filename1, const std::string& filename2){
    this->deleteWriter();
    this->writer1 = new Writer(filename1, this->opt->compression);
    this->writer2 = new Writer(filename2, this->opt->compression);
}

void ThreadConfig::initWriter(std::ofstream* stream){
    this->deleteWriter();
    this->writer1 = new Writer(stream);
}

void ThreadConfig::initWriter(std::ofstream* stream1, std::ofstream* stream2){
    this->deleteWriter();
    this->writer1 = new Writer(stream1);
    this->writer2 = new Writer(stream2);
}

void ThreadConfig::initWriter(gzFile gzfile){
    this->deleteWriter();
    this->writer1 = new Writer(gzfile);
}

void ThreadConfig::initWriter(gzFile gzfile1, gzFile gzfile2){
    this->deleteWriter();
    this->writer1 = new Writer(gzfile1);
    this->writer2 = new Writer(gzfile2);
}

void ThreadConfig::addFilterResult(int result){
    this->filterResult->addFilterResult(result);
}

void ThreadConfig::initWriterForSplit(){
    if(this->opt->out1.empty()){
        return;
    }
    std::string num = std::to_string(this->workingSplit + 1);
    if(this->opt->digits > 0){
        while(num.size() < this->opt->digits){
            num = "0" + num;
        }
    }

    std::string filename1 = util::joinpath(util::dirname(this->opt->out1), num + "." + util::basename(this->opt->out1));
    if(this->opt->paired){
        this->initWriter(filename1);
    }else{
        std::string filename2 = util::joinpath(util::dirname(this->opt->out2), num + "." + util::basename(this->opt->out2));
    }
}

void ThreadConfig::markProcessed(size_t readNum){
    this->currentSplitReads += readNum;
    if(!this->opt->enableSplit){
        return;
    }
    if(this->currentSplitReads >= this->opt->splitSize){

    }

#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "util.h"
#include "stats.h"
#include "writer.h"
#include "filterresult.h"

struct ThreadOpt{
    int estimatedReadLen;
    bool enableSplit;
    bool splitByFileNumber;
    bool splitByFileLines;
    int thread;
    int splitNumber;
    size_t splitSize;
    bool paired;
    int compression;
    int digits;
    std::string out1;
    std::string out2;
};

class ThreadConfig{
    public:
        ThreadOpt* opt;
        FilterOpt* fopt;
        Stats* preStats1;
        Stats* preStats2;
        Stats* postStats1;
        Stats* postStats2;
        Writer* writer1;
        Writer* writer2;
        FilterResult* filterResult;

        int threadId;
        int workingSplit;
        size_t currentSplitReads;
        bool canbeStopped;

    public:
        ThreadConfig(ThreadOpt* opt, int threadId, bool paired = false);
        ~ThreadConfig();

        inline Stats* getPreStats1(){return this->preStats1;}
        inline Stats* getPreStats2(){return this->preStats2;}
        inline Stats* getPostStats1(){return this->postStats1;}
        inline Stats* getPostStats2(){return this->postStats2;}
        inline Writer* getWriter1(){return this->writer1;}
        inline Writer* getWriter2(){return this->writer2;}
        inline FilterResult* getFilterResult(){return this->filterResult;}

        void initWriter(const std::string& filename);
        void initWriter(const std::string& filename1, const std::string& filename2);
        void initWriter(std::ofstream* stream);
        void initWriter(std::ofstream* stream1, std::ofstream* stream2);
        void initWriter(gzFile gzfile);
        void initWriter(gzFile gzfile1, gzFile gzfile2);

        void addFilterResult(int result);
        inline int getThreadId(){return this->threadId;}
        
        void markProcessed(size_t readNum);
        void initWriterForSplit();
        bool canBeStopped();
        void cleanup();

        void deleteWriter();
        void writeEmptyFilesForSplitting();
};
        
#endif

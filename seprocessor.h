#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <memory>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <atomic>
#include <unistd.h>
#include "read.h"
#include "filter.h"
#include "fqreader.h"
#include "util.h"
#include "adaptertrimmer.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "polyx.h"
#include "threadconfig.h"
#include "duplicate.h"
#include "stats.h"
#include "filterresult.h"

namespace fqlib{
    struct ReadPack{
        Read** data;
        int count;
    };

    struct ReadRepository{
        ReadPack** packBuffer;
        std::atomic<long> readPos;
        std::atomic<long> writePos;
    };

    class SingleEndProcessor{
        public:
            SingleEndProcessor(Options* opt);
            ~SingleEndProcessor();
            bool process();
            
            bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
            void initPackRepository();
            void destroyPackRepository();
            void producePack(ReadPack* pack);
            void consumePack(ReadPack* pack);
            void producerTask();
            void consumerTask(ThreadConfig* config);
            void initConfig(ThreadConfig* config);
            void initOutput();
            void closeOutput();
            void writeTask(WriterThread* config);
            
            Options* mOptions;
            ReadRepository mRepo;
            std::atomic<bool> mProduceFinished;
            std::atomic<int> mFinishedThreads;
            std::mutex mInputMtx;
            std::mutex mOutputMtx;
            Filter* mFilter;
            gzFile mZipFile;
            std::ofstream* mOutStream;
            UmiProcessor* mUmiProcessor;
            WriterThread* mLeftWriter;
            Duplicate* mDuplicate;
    };
}

#endif

#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <mutex>
#include <atomic>
#include "common.h"
#include "util.h"
#include "writer.h"
#include "options.h"

namespace fqlib{
    class WriterThread{
        public:
            WriterThread(Options* opt, const std::string& filename);
            ~WriterThread();
    
            void iniWriter(const std::string& filename);
            void iniWriter(std::ofstream* ofs);
            void iniWriter(gzFile gzfile);
    
            void cleanup();
            bool isCompleted();
            void output();
            void input(char* cstr, size_t size);
            bool setInputCompleted();
            
            size_t bufferLength();
            
            /** get output fiilename
             * @return filename
             */
            inline std::string getFilename(){
                return mFilename;
            }
    
            void deleteWriter();
            
            Options* mOptions;
            Writer* mWriter;
            std::string mFilename;
    
            bool mInputCompleted;
            std::atomic_long mInputCounter;
            std::atomic_long mOutputCounter;
            char** mRingBuffer;
            size_t* mRingBufferSizes;
    
            std::mutex mtx;
    };
}
#endif

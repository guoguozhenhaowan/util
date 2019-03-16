#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <mutex>
#include <atomic>
#include "util.h"
#include "writer.h"

class WriterThread{
    public:
        WriterThread(const std::string& filename, int compression);
        ~WriterThread();

        void iniWriter(const std::string& filename, int compression);
        void iniWriter(std::ofstream* ofs);
        void iniWriter(gzFile gzfile);

        void cleanup();
        bool isCompleted();
        void output();
        void input(char* cstr, size_t size);
        bool setInputCompleted();
        
        size_t bufferLength();
        
        /** get output filename
         * @return outpuf filename
         */
        inline std::string getFilename(){
            return this->filename;
        }

        void deleteWriter();

        Writer* writer;
        std::string filename;

        bool inputCompleted;
        std::atomic_long inputCounter;
        std::atomic_long outputCounter;
        char** ringBuffer;
        size_t* ringBufferSizes;

        std::mutex mtx;
};

#endif

#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <cstdio>
#include <cstlib>
#include <string>
#include <memory>
#include <cmath>
#include "overlapanalysis.h"
#include "options.h"
#include "read.h"
#include "common.h"

namespace fqlib{
    class Duplicate{
        public:
            Duplicate(Options* opt);
            ~Duplicate();
            
            void statRead(Read* r1);
            void statPair(Read* r1, Read* r2);
            void addRecord(int32_t key, u_int64_t kmer32, u_int8_t gc);
            double statAll(int* hist, double* meanGC, int histSize);
        
        private:
            Options* mOptions;
            int mKeyLenInBase;
            int mKeyLenInBit;
            u_int64_t* mDups;
            u_int64_t* mCounts;
            u_int8_t* mGC;
    };
}
        
#endif

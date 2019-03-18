#ifndef POLY_X_H
#define POLY_X_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "filterresult.h"

namespace fqlib{
    /** Class to do Poly X trimming */
    class PolyX{
        public:
            /** construct a PolyX object */
            PolyX();
            /** Destroy a PolyX object */
            ~PolyX();

            /** trim polyG from 3' end\n
             * if at most 5 mismatch in length l, and l >= compareReq, trim it
             * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
             * @param r1 pointer to Read object(read1)
             * @param r2 pointer to Read object(read2)
             * @param fr pointer to FilterResult object(may be NULL)
             * @param compareReq required length of sequence to be polyG
             */
            static void trimPolyG(Read* r1, Read* r2, FilterResult* fr, int compareReq);
            
            /** trim polyG from 3' end\n
             * if at most 5 mismatch in length l, and l >= compareReq, trim it
             * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
             * @param r pointer to Read object
             * @param fr pointer to FilterResult object(may be NULL)
             * @param compareReq required length of sequence to be polyG
             */
            static void trimPolyG(Read* r, FilterResult* fr, int compareReq);
            
            /** trim polyX from 3' end\n
             * if at most 5 mismatch in length l, and l >= compareReq, trim it
             * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
             * @param r1 pointer to Read object(read1)
             * @param r2 pointer to Read object(read2)
             * @param fr pointer to FilterResult object(may be NULL)
             * @param compareReq required length of sequence to be polyG
             */
            static void trimPolyX(Read* r1, Read* r2, FilterResult* fr, int compareReq);

            /** trim polyX from 3' end\n
             * if at most 5 mismatch in length l, and l >= compareReq, trim it
             * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
             * @param r pointer to Read object
             * @param fr pointer to FilterResult object(may be NULL)
             * @param compareReq required length of sequence to be polyG
             */
            static void trimPolyX(Read* r, FilterResult* fr, int compareReq);
    };
}

#endif

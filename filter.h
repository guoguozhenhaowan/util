#ifndef FILTER_H
#define FILTER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "common.h"
#include "read.h"


/** Struct to store various quality threshold dependent fastq read cut/trim options */
struct QualCutOpt{
    int minFrontQual;  ///< minimum average window quality required to stop sliding cut from front(5'end)
    int minTailQual;   ///< minimum average window quality required to stop sliding cut from tail(3'end)
    int minRightQual;  ///< minimum average window quality required to stop sliding detect a low quality window from front(5'end)
    int windowFront;   ///< window size for sliding cut from front 
    int windowTail;    ///< window size for sliding cut from tail
    int windowRight;   ///< window size for sliding detect low quality window from front
    bool cutFront;     ///< if true, sliding window from front(5'end) until average quality > minFrontQual and cut from the last non-N base of this window
    bool cutTail;      ///< if true, sliding window from tail(3'end) until average quality > minTailQual and cut from the last non-N base of this window
    bool cutRright;    ///< if true, sliding window from front(5'end) until average quality < minRightQual and cut before the first base with quality < minRightQual
                       ///< this option is excluded with cutTail
};

/** Class to do fastq read filter by various standards and methods */
class Filter{
        int lowQualThreshold;    ///< low quality threshold, a base with quality greater than lowQualThreshold is not low quality base
        int maxLowQualNum;       ///< maximum low quality bases allowed for a sequence
        int maxNBaseNum;         ///< maximum N bases allowed for a sequence
        double minComplexity;    ///< minimum complexity allowed for a sequence
        int minLen;              ///< minimum length required for a sequence
        int maxLen;              ///< maximum length allowed for a sequence
    public:
        /** Construct a Filter object, negative parameter will turn the corresponding filterr
         * @param lowQualThreshold low quality threshold, if set to negative number
         * @param maxLowQualNum maximum low quality bases allowed for a sequence
         * @param maxNBaseNum maximum N bases allowed for a sequence
         * @param minComplexity minimum complexity allowed for a sequence
         * @param minLen minimum length required for a sequence
         * @param maxLen maximum length allowed for a sequence
         */
        Filter(int lowQualThreshold, int maxLowQualNum, int maxNBaseNum, double minComplexity, int minLen, int maxLen){
            this->lowQualThreshold = lowQualThreshold;
            this->maxLowQualNum = maxLowQualNum;
            this->maxNBaseNum = maxNBaseNum;
            this->minComplexity = minComplexity;
            this->minLen = minLen;
            this->maxLen = maxLen;
        }

        /** Destroy a Filter object */
        ~Filter() = default;

        /** Test whether a Read pass filter
         * @param r pointer to a Read object
         * @return 0 if passed
         */
        int passFilter(Read* r);
        /** Test whether a Read pass low complexity filter
         * @param r pointer to a Read object
         * @return true if r is not low complexity
         */
        bool passLowComplexityFliter(Read* r);

        /** Trim and cut a Read by various methods
         * @param r pointer to a Read object
         * @param forceFrontCut force cut length from 5'end, will always happen
         * @param forceTailCut force cut length from 3'end, will always happen
         * @param opt vasious quality cut options, some may happen
         * @return a Read trim and cut properly
         */
        static Read* trimAndCut(Read* r, int forceFrontCut, int forceTailCut, QualCutOpt* opt);
        
        /** Filter a Read by its first index
         * @param r pointer to a Read object
         * @param indexBlackList a list of index to search against
         * @param threshold maximum different bases allowed for a succesful match
         * @return true if r's index is succesfully matched in the indexBlackList
         */
        inline static bool filterByIndex(Read* r, const std::vector<std::string>& indexBlackList, const int& threshold){
            return Filter::match(indexBlackList, r->firstIndex(), threshold);
        }
        
        /** Filter a pair of Read by their first index
         * @param r1 pointer to a Read object
         * @param r2 pointer to a Read object
         * @param indexBlackList a list of index to search against
         * @param threshold maximum different bases allowed for a succesful match
         * @return true if the first index of r1 or r2 succesfully matched in the indexBlackList
         */
        inline static bool filterByIndex(Read* r1, Read* r2, const std::vector<std::string>& indexBlackList, const int& threshold){
            return Filter::filterByIndex(r1, indexBlackList, threshold) || Filter::filterByIndex(r2, indexBlackList, threshold);
        }

        /** match a sequnce (target) against a series of sequences (list)
         * @param list a vector of strings to be matched against
         * @param target a string to be tested whether it will match any string in list 
         * @param threshold maximum character difference allowed for a successful match
         * @return true if any successful match occurs between target and string in list
         */
        static bool match(const std::vector<std::string>& list, const std::string& target, int threshold);
};

#endif

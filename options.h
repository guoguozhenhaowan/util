#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

/** struct to store various quality threshold dependent fastq read cut/trim options */
struct QualityCutOptions{
    bool enableFront;      ///< if true, sliding window from front(5'end) until average quality > minFrontQual and cut from the last non-N base of this window
    bool enableTail;       ///< if true, sliding window from tail(3'end) until average quality > minTailQual and cut from the last non-N base of this window
    bool enableRright;     ///< if true, sliding window from front(5'end) until average quality < minRightQual and cut before the first base with quality < minRightQual
    int qualityShared;     ///< shared low quality threshold for cutting, make setting default values easy
    int windowSizeShared;  ///< shared window size for sliding, make setting default values easy
    int qualityFront;      ///< minimum average window quality required to stop sliding cut from front(5'end)
    int qualityTail;       ///< minimum average window quality required to stop sliding cut from tail(3'end)
    int qualityRight;      ///< minimum average window quality required to stop sliding detect a low quality window from front(5'end)
    int windowSizeFront;   ///< window size for sliding cut from front
    int windowSizeTail;    ///< window size for sliding cut from tail
    int windowSizeRight;   ///< window size for sliding detect low quality window from front
    
    /** construct a QualityCutOptions object and set default values */
    QualityCutOptions(){
        enableFront = false; 
        enableTail = false;  
        enableRright = false;
        qualityShared = 20;       
        windowSizeShared = 4;    
        qualityFront = qualityShared;     
        qualityTail = qualityShared;         
        qualityRight = qualityShared;        
        windowSizeFront = windowSizeShared;     
        windowSizeTail = windowSizeShared;      
        windowSizeRight = windowSizeShared;
    }
};

/** struct to store index filtering options */
struct IndexFilterOptions{
    bool enabled;                         ///< enable index filtering if set true
    int threshold;                        ///< maximum different bases allowed for an external index hit blacklist successfully 
    std::vector<std::string> blacklist1;  ///< read1 index blacklist
    std::vector<std::string> blacklist2;  ///< read2 index blacklist

    /** construct a IndexFilterOptions object and set default values */
    IndexFilterOptions(){
        enabled = false;
        threshold = 0;
    }
};

/** struct to store over representation sequence analyisis options */
struct OverrepresentedSequenceAnalysisOptions{
    bool enabled;   ///< enable over representation sequence analyisis
    int sampling;   ///< sampling frequence for ORA
    
    /** construct a OverrepresentedSequenceAnalysisOptions object and set default values */
    OverrepresentedSequenceAnalysisOptions(){
        enabled = false;
        sampling = 20;
    }
};

/** struct to store base correction options */
struct  CorrectionOptions {
    bool enabled; ///< enable base correction if true
    /** construct a CorrectionOptions object and set default values */
    CorrectionOptions() {
        enabled = false;
    }
};

/** struct to store low complexity filter options */
struct LowComplexityFilterOptions {
    bool enabled;     ///< enable low complexity filter if true
    double threshold; ///< threshold to test sequence complexity, below which will be dropped 
    /** construct a LowComplexityFilterOptions object and set default values */
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
};

/** struct to store read length filter options */
struct ReadLengthFilterOptions{
    bool enabled;         ///< enable read length filter if true
    int minReadLength;    ///< if read_length < minReadLength, then this read will be discarded
    int maxReadLength;    ///< if read_length > maxReadLength, then this read will be discarded
    /** construct a ReadLengthFilterOptions object and set default values */
    ReadLengthFilterOptions(){
        enabled = false;
        minReadLength = 15;
        // 0 for no limination
        maxReadLength = 0; 
    }
}

/** struct to store adapter detect/trimming options */
struct AdapterOptions{
    bool enableTriming;               ///< enable index trimming if true
    bool enableDetectForPE;           ///< enable auto detection of index for pair end reads if true
    bool adapterSeqR1Provided;        ///< adapter sequence for read1 is provided externally if true
    bool adapterSeqR2Provided;        ///< adapter sequence for read2 is provided externally if true
    std::string inputAdapterSeqR1;    ///< adapter sequence for read1 provided externally
    std::string inputAdapterSeqR2;    ///< adapter sequence for read2 provided externally
    std::string detectedAdapterSeqR1; ///< adapter sequence for read1 auto detected 
    std::string detectedAdapterSeqR2; ///< adapter sequence for read2 auto detected
    /** construct a AdapterOptions object and set default values */
    AdapterOptions(){
        enableTriming = true;
        adapterSeqR1Provided = false;
        adapterSeqR2Provided = false;
        enableDetectForPE = true;
    }
};

#endif

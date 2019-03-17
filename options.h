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
    bool enabled;                              ///< enable over representation sequence analyisis
    int sampling;                              ///< sampling frequence for ORA
    std::map<std::string, long> overRepSeqR1;  ///< over represented sequences count of read1
    std::map<std::string, long> overRepSeqR2;  ///< over represented sequences count of read2
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
};

/** struct to store quality filter options */
struct QualityFilterOptions{
    bool enabled;               ///< quality filter enabled if true
    char lowQualityLimit;       ///< if a base's quality < lowQualityLimit, hen it's considered as a low quality base
    int lowQualityBaseLimit;    ///< if low quality bases number > lowQualityBaseLimit, then discard this read
    int nBaseLimit;             ///< if N bases number > nBaseLimit, then discard this read
    /** construct a QualityFilterOptions object and set default values */
    QualityFilterOptions(){
        enabled = true;
        // '0' == Q15
        lowQualityLimit = '0';
        lowQualityBaseLimit = 40;
        nBaseLimit = 5;
    }
};

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

/** struct to store force trimming options */
struct ForceTrimOptions{
    int front1;  ///< first cycles trimmed for read1
    int tail1;   ///< last cycles trimmed for read1
    int front2;  ///< first cycles trimmed for read2
    int tail2;   ///< last cycles trimmed for read2
    int maxLen1; ///< maximum length of read1
    int maxLen2; ///< maximum length of read2
    /** construct a ForceTrimOptions object and set default values */
    ForceTrimOptions(){
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
};

/** struct to hold various option structs and interface  options together */
struct Options{
    // program interface options 
    std::string in1;        ///< input read1 filename
    std::string in2;        ///< input read2 filename
    std::string out1;       ///< output read1 filename
    std::string out2;       ///< output read2 filename
    std::string jsonFile;   ///< output json filename
    std::string htmlFile;   ///< output html report filename
    std::string reportTitle; ///< html report title
    int compression;         ///< compression level for gz format output
    bool phred64;            ///< the input file is using phred64 quality scoring if true 
    bool donotOverwrite;     ///< do not over write existing files
    bool inputFromSTDIN;     ///< read from STDIN
    bool outputToSTDOUT;     ///< write to STDOUT
    bool interleavedInput;   ///< the input read1(in1) file is an interleaved PE fastq
    int readsToProcess;      ///< number of reads to process
    int thread;              ///< number of threads to do paralel work
    bool verbose;            ///< output debug information if true
    // common options
    int seqLen1;             ///< estimated maximum read length of read1
    int seqLen2;             ///< estimated maximum read length of read2
    int insertSizeMax;       ///< maximum value of insert size
    int overlapRequire;      ///< overlap region minimum length
    int overlapDiffLimit;    ///< overlap region maximum different bases allowed
    // submodule options
    ForceTrimOptions trim;                ///< ForceTrimOptions object
    QualityFilterOptions qualFilter;      ///< QualityFilterOptions object
    QualityCutOptions qualitycut;         ///< QualityCutOptions object 
    ReadLengthFilterOptions lengthFilter; ///< ReadLengthFilterOptions object 
    AdapterOptions adapter;               ///< AdapterOptions object 
    CorrectionOptions correction;         ///< CorrectionOptions object 
    OverrepresentedSequenceAnalysisOptions overRepAna; ///< OverrepresentedSequenceAnalysisOptions object
    LowComplexityFilterOptions complexityFilter;       ///< LowComplexityFilterOptions object
    IndexFilterOptions indexFilter;                    ///< IndexFilterOptions object
};

#endif

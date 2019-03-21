#include "peprocessor.h"

namespace fqlib{
    PairEndProcessor::PairEndProcessor(Options* opt){
        mOptions = opt;
        mProduceFinished = false;
        mFinishedThreads = 0;
        mFilter = new Filter(opt);
        mOutStream1 = NULL;
        mOutStream2 = NULL;
        mZipFile1 = NULL;
        mZipFile2 = NULL;
        mUmiProcessor = new UmiProcessor(mOptions);

        int insertSizeBufLen = mOptions->insertSizeMax + 1;
        mInsertSizeHist = new long[insertSizeBuflen];
        std::memset(mInsertSizeHist, 0, sizeof(long) * insertSizeBufLen);
        mLeftWriter = NULL;
        mRightWriter = NULL;

        mDuplicate = NULL;
        if(mOptions->duplicate.enabled){
            mDuplicate = new Duplicate(mOptions);
        }
    }

    PairEndProcessor::~PairEndProcessor(){
        delete mInsertSizeHist;
        delete mUmiProcessor;
        if(mDuplicate){
            delete mDuplicate;
            mDuplicate = NULL;
        }
    }

    void PairEndProcessor::initOutput(){
        if(mOptions->out1.empty() || mOptions->out2.empty()){
            return;
        }
        mLeftWriter = new WriterThread(mOptions, mOptions->out1);
        mRightWriter = new WriterThread(mOptions, mOptions->out2);
    }

    void PairEndProcessor::closeOutput(){
        if(mLeftWriter){
            delete mLeftWriter;
            mLeftWriter = NULL;
        }
        if(mRightWriter){
            delete mRightWriter;
            mRightWriter = NULL;
        }
    }

    void PairEndProcessor::initConfig(ThreadConfig* config){
        if(mOptions->out1.empty()){
            return;
        }
        if(mOptions->split.enabled){
            config->initWriterForSplit();
        }
    }

    bool PairEndProcessor::process(){
        if(!mOptions->split.enabled){
            initOutput();
        }

        initReadPairPackRepository();
        std::thread producer(


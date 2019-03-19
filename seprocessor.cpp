#include "seprocessor.h"

namespace fqlib{
    SingleEndProcessor::SingleEndProcessor(Options* opt){
        mOptions = opt;
        mProduceFinished = false;
        mFinishedThreads = 0;
        mFilter = new Filter(mOptions);
        mOutStream = NULL;
        mZipFile = NULL;
        mUmiProcessor = new UmiProcessor(mOptions);
        mLeftWriter = NULL;
        mDuplicate = NULL;
        if(mOptions->duplicate.enabled){
            mDuplicate = new Duplicate(mOptions);
        }
    }

    SingleEndProcessor::~SingleEndProcessor(){
        delete mFilter;
        delete mUmiProcessor;
        if(mDuplicate){
            delete mDuplicate;
            mDuplicate = NULL;
        }
    }

    void SingleEndProcessor::initOutput(){
        if(mOptions->out1.empty()){
            return;
        }
        mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    }

    void SingleEndProcessor::closeOutput(){
        if(mLeftWriter){
            delete mLeftWriter;
            mLeftWriter = NULL;
        }
    }

    void SingleEndProcessor::initConfig(ThreadConfig* config){
        if(mOptions->out1.empty()){
            return;
        }
        
        if(mOptions->split.enabled){
            config->initWriterForSplit();
        }
    }

    void SingleEndProcessor::initPackRepository(){
        mRepo.packBuffer = new ReadPack*[compar::PACK_NUM_LIMIT];
        std::memset(mRepo.packBuffer, 0, sizeof(ReadPack*) * compar::PACK_NUM_LIMIT);
        mRepo.writePos = 0;
        mRepo.readPos = 0;
    }

    void SingleEndProcessor::producePack(ReadPack* pack){
        mRepo.packBuffer[mRepo.writePos] = pack;
        ++mRepo.writePos;
    }

    void SingleEndProcessor::producerTask(){
        if(mOptions->verbose){
            util::loginfo("start to load data");
        }
        long lastReported = 0; // total number of reads have been loaded into memory
        int slept = 0;
        long readNum = 0;  // total number of reads have been loaded into memory and put into mReop
        bool splitSizeReEvaluated = false;
        Read** data = new Read*[compar::PACK_SIZE];
        std::memset(data, 0, sizeof(Read*) * compar::PACK_SIZE);
        FqReader reader(mOptions->in1, true, mOptions->phred64);
        int count = 0; // number of reads have been loaded into memory and put into data but not mReop 
        bool needToBreak = false;
        while(true){
            Read* read = reader.read();
            if(!read || needToBreak){
                // last read have got or the first N reads needed have got
                ReadPack* pack = new ReadPack;
                pack->data = data;
                pack->count = count;
                producePack(pack);
                data = NULL;
                if(read){
                    delte read;
                    read = NULL;
                }
                break;
            }
            data[count] = read;
            ++count;
            // if only need to process the first N reads
            if(mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess){
                needToBreak = true;
            }
            // if verbose log #M reads loaded into memory
            if(mOptions->verbose && count + readNum >= lastReported + 1000000){
                lastReported = count + readNum;
                std::string msg = "loaded " + std::to_string(lastReported/1000000) + "M reads";
                util::loginfo(msg);
            }
            // the pack is full or first N reads needed have got
            if(count == compar::PACK_SIZE || needToBreak){
                ReadPack* pack = new ReadPack;
                pack->data = data;
                pack->count = count;
                producePack(pack);
                //re-initialize data for next pack
                data = new Read*[compar::PACK_SIZE];
                std::memset(data, 0, sizeof(Read*) * compar::PACK_SIZE);
                // if the consumer is far behind this producer, sleep and wait to limit memory usage
                while(mRepo.writePos - mRepo.readPos > compar::PACK_IN_MEM_LIMIT){
                    ++slept;
                    usleep(100);
                }
                readNum += count;
                // if the writer threads are far behind this producer, sleep and wait
                if(readNum % (compar::PACK_SIZE * compar::PACK_IN_MEM_LIMIT) == 0 && mLeftWriter){
                    while(mLeftWriter->bufferLength() > compar::PACK_IN_MEM_LIMIT){
                        ++slept;
                        usleep(100);
                    }
                }
                // rest count to 0
                count = 0;
            }
        }
        mProduceFinished = true;
        if(mOptions->verbose){
            std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " finished";
            util::loginfo(msg);
        }
        if(data != NULL){
            delete[] data;
        }
    }
    
    void SingleEndProcessor::consumerTask(ThreadConfig* config){
        while(true){
            if(config->canBeStopped){
                ++mFinishedThreads;
                break;
            }
            while(mRepo.writePos <= mRepo.readPos){
                if(mProduceFinished){
                    break;
                }
                usleep(1000);
            }
            if(mProduceFinished && mRepo.writePos == mRepo.readPos){
                ++mFinishedThreads;
                if(mOptions->verbose){
                    std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " data processing finished";
                    util::loginfo(msg);
                }
                break;
            }
            if(mProduceFinished){
                if(mOptions->verbose){
                    std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " is processing the " +
                        std::string(mRepo.readPos) + "/" + std::to_string(mRepo.writePos) + " pack";
                    util::loginfo(msg);
                }
                consumePack(config);
            }
        }

        if(mFinishedThreads == mOptions->thread){
            if(mLeftWriter){
                mLeftWriter->setInputCompleted();
            }
        }

        if(mOptions->verbose){
            std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + 
                " finished";
            util::loginfo(msg);
        }
    }

    void SingleEndProcessor::consumePack(ThreadConfig* config){
        ReadPack* data;
        mInputMtx.lock();
        while(mRepo.writePos <= mRepo.readPos){
            usleep(1000);
            if(mProduceFinished){
                mInputMtx.unlock();
                return;
            }
        }
        data = mRepo.packBuffer[mRepo.readPos];
        ++mRepo.readPos;
        mInputMtx.unlock();
        processSingleEnd(data, config);
    }


    bool SingleEndProcessor::process(){
        if(!mOptions->split.enabled){
            initOutput();
        }

        initPackRepository();
        std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

        int cycle = 151;
        ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
        for(int t = 0; t < mOptions->thread; ++t){
            configs[t] = new ThreadConfig(mOptions, t, false);
            initConfig(configs[t]);
        }

        std::thread** threads = new std::thread*[mOptions->thread];
        for(int t = 0; t < mOptions->thread; ++t){
            threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
        }

    }

}

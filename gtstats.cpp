#include "stats.h"
#include "gtest/gtest.h"

class StatsTest : public testing::Test{
    protected:
        void SetUp() override{
            s.setKmerLen(0);
            s.setOverRepSampleFreq(0);
            s.allocateRes();
        }
        Evaluator e = {"./testdata/R1.adaptor.fq.gz", 0};
        FqReader f = {"./testdata/R1.adaptor.fq.gz"};
        Stats s = {e.getReadLen()};
};

TEST_F(StatsTest, statRead){
    const int maxReadInMemory = 10000;
    std::vector<Stats*> l;
    l.reserve(maxReadInMemory + 1);
    l.resize(1);
    Read* r = NULL;
    Stats* fs = NULL;
    int readInMemory = 0;
    while((r = f.read()) != NULL){
        Stats* ts = new Stats(e.getReadLen());
        ts->setKmerLen(6);
        ts->setOverRepSampleFreq(100);
        ts->allocateRes();
        ts->statRead(r);
        ts->summarize();
        l.push_back(ts);
        ++readInMemory;
        delete r;
        if(readInMemory >= maxReadInMemory){
            if(!fs){
                l[0] = l[l.size()-1];
                l.resize(l.size() -1);
            }else{
                l[0] = fs;
            }
            fs = Stats::merge(l);
            for(size_t i = 1; i < l.size(); ++i){
                delete l[i];
            }
            l.resize(1);
            break;
            readInMemory = 0;
        }
    }
    std::ofstream fw("test.json");
    fs->reportJson(fw, " ");
    fw.close();
    fw.open("test.html");
    fs->reportHtml(fw, "Before Filtering", "TestLibrary");
    fw.close();
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

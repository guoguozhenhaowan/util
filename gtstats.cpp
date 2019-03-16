#include "stats.h"
#include "gtest/gtest.h"

class StatsTest : public testing::Test{
    protected:
        void SetUp() override{
            s.setKmerLen(6);
            s.setOverRepSampleFreq(100);
            s.allocateRes();
            s.initOverRepSeq(e); // must do initOverRepSeq before doing ORA
        }
        Evaluator e = {"./testdata/R1.adaptor.fq.gz", 0};
        FqReader f = {"./testdata/R1.adaptor.fq.gz"};
        Stats s = {e.getReadLen()};
};

TEST_F(StatsTest, statRead){
    Read* r = NULL;
    while((r = f.read()) != NULL){
        s.statRead(r);
        delete r;
    }
    s.summarize();
    std::ofstream fw("test.json");
    s.reportJson(fw, " ");
    fw.close();
    fw.open("test.html");
    s.reportHtml(fw, "Before Filtering", "TestLibrary");
    fw.close();
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

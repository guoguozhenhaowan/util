#include "stats.h"
#include "gtest/gtest.h"

TEST(StatsTest, statRead){
    Options opt;
    opt.in1 = "./testdata/R1.adaptor.fq.gz";
    opt.overRepAna.enabled = true;
    opt.overRepAna.sampling = 100;

    Evaluator e(&opt);
    e.computeOverRepSeq(opt.in1, opt.overRepAna.overRepSeqR1);
    e.evaluateReadLen();
    opt.kmer.enabled = true;
    opt.kmer.kmerLen = 4;
    Stats s(&opt, false);
    
    FqReader f = {"./testdata/R1.adaptor.fq.gz"};
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

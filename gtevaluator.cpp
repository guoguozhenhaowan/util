#include "gtest/gtest.h"
#include "evaluator.h"

TEST(EvaluatorTest, All){
    fqlib::Options opt = {};
    opt.in1 = "./testdata/R1.adaptor.fq.gz";
    opt.trim.tail1 = 0;
    fqlib::Evaluator e = {&opt};
    e.evaluateReadLen();
    e.evaluateReadNum();
    e.evaluateAdapterSeq(false);
    e.evaluateTwoColorSystem();
    
    EXPECT_EQ(opt.est.adapter, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    EXPECT_TRUE(opt.est.illuminaAdapter);
    EXPECT_TRUE(opt.est.twoColorSystem);
    EXPECT_EQ(opt.est.seqLen1, 150);
    EXPECT_TRUE(std::abs((long)opt.est.readsNum - 300000) < 1000);
    std::cout << "estimated readnum: " << opt.est.readsNum << std::endl;
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include "filter.h"
#include "gtest/gtest.h"

TEST(Filter, trimAndCut){
    Read r("@name", 
           "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT", 
           "+",
           "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    QualCutOpt opt;
    opt.cutFront = true;
    opt.cutTail = true;
    opt.windowFront = 4;
    opt.minFrontQual = 20;
    opt.windowTail = 4;
    opt.minTailQual = 20;

    Read* ret = Filter::trimAndCut(&r, 0, 1, &opt);

    EXPECT_EQ(ret->seq.seqStr, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    EXPECT_EQ(ret->quality, "CCCCCCCCCCC////CCCCCCCCCCCCC");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

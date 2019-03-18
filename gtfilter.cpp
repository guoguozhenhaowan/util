#include "filter.h"
#include "gtest/gtest.h"

TEST(Filter, trimAndCut){
    fqlib::Read r("@name", 
           "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT", 
           "+",
           "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    fqlib::Options opt;
    opt.qualitycut.enableFront = true;
    opt.qualitycut.enableTail = true;
    opt.qualitycut.windowSizeFront = 4;
    opt.qualitycut.qualityFront = 20;
    opt.qualitycut.windowSizeTail = 4;
    opt.qualitycut.qualityTail = 20;

    fqlib::Filter f(&opt);
    fqlib::Read* ret = f.trimAndCut(&r, 0, 1);

    EXPECT_EQ(ret->seq.seqStr, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    EXPECT_EQ(ret->quality, "CCCCCCCCCCC////CCCCCCCCCCCCC");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

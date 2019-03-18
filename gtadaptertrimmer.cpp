#include "gtest/gtest.h"
#include "adaptertrimmer.h"

TEST(AdapterTrimmer, trimBySequence){
    fqlib::Read r("@name",
                  "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
                  "+",
                  "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    fqlib::Read r2 = r;

    std::string adapter = "TTTTCCACGGGGATACTACTG";
    bool trimmed = fqlib::AdapterTrimmer::trimBySequence(&r, NULL, adapter);
    EXPECT_EQ(r.seq.seqStr, "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");
    
    adapter = "CATTATTTAACCCCCCCCCCCCCCCCCCCCCCCC";
    trimmed = fqlib::AdapterTrimmer::trimBySequence(&r2, NULL, adapter);
    EXPECT_EQ(r2.seq.seqStr, "");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

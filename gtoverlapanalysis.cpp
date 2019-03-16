#include "gtest/gtest.h"
#include "overlapanalysis.h"

TEST(OverlapAnalysis, analyze){
    // TEMP_LEN > SEQ_LEN
    Seq s1("CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC");
    Seq s2("ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC");
    OverlapResult ovr = ::OverlapAnalysis::analyze(s1, s2);
    EXPECT_EQ(ovr.overlapLen, 79);
    EXPECT_TRUE(ovr.overlapped);
    EXPECT_EQ(ovr.offset, 10);
    EXPECT_EQ(ovr.diff, 1);
    
    // TEMP_LEN < SEQ_LEN
    Seq sa("NGAATCTTTCTCTGGAAGAGGCTCATATGGTGGGGCAAATACTGGACCTTTATGTTCTAGGAATTTCCACTTGATGCCTTCAGGATAGCGCTCTTCTTCCCACCATTTCCACTTCTGTTCCTCTTCTTTAAACAGGCGAGATCGGAAGAG");
    Seq sb("CGCCTGTTTAAAGAAGAGGAACAGAAGTGGAAATGGTGGGAAGAAGAGCGCTATCCTGAAGGCATCAAGTGGAAATTCCTAGAACATAAAGGTCCAGTATTTGCCCCACCATATGAGCCTCTTCCAGAGAAAGATTCCAGATCGGAAGAG");

    ovr = OverlapAnalysis::analyze(sa, sb);
    EXPECT_EQ(ovr.overlapLen, 138);
    EXPECT_TRUE(ovr.overlapped);
    EXPECT_EQ(ovr.offset, -12);
    EXPECT_EQ(ovr.diff, 1);
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
#include "read.h"
#include "gtest/gtest.h"

class ReadTest : public ::testing::Test{
    protected:
        void SetUp() override{
        }

        fqlib::Read r = {"@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
               "CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC",
               "+",
               "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA"};
        fqlib::Read r_one_index = {"@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT",
               "CTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTCCTTAGGAGGACATTTTTTACATGAAATTATTAACCTAAATAGAGTTGATC",
               "+",
               "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEEEAEEEEEEEEEEEEEEEAEEE/EEEEEEEEEEAAEAEAAEEEAEEAA"};
};

TEST_F(ReadTest, reverseComplement){
    fqlib::Read* pR = r.reverseComplement();
    EXPECT_EQ(pR->seq.seqStr, r.seq.reverseComplement().seqStr); 
    EXPECT_EQ(pR->strand, "-");
    EXPECT_EQ(pR->name, r.name);
    EXPECT_EQ(pR->quality, std::string(r.quality.rbegin(), r.quality.rend()));
}

TEST_F(ReadTest, firstIndex){
    EXPECT_EQ(r.firstIndex(), "TATAGCCT");
    EXPECT_EQ(r_one_index.firstIndex(), "TATAGCCT");
}

TEST_F(ReadTest, lastIndex){
    EXPECT_EQ(r.lastIndex(), "GGTCCCGA");
    EXPECT_EQ(r_one_index.lastIndex(), "TATAGCCT");
}

TEST_F(ReadTest, length){
    EXPECT_EQ(r.length(), 151);
}

class ReadPairTest : public ::testing::Test{
    protected:
        void SetUp() override{
            rp.left = new fqlib::Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                               "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAG",
                               "+",
                               "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
            rp.right = new fqlib::Read("@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA",
                                "AAAAAACTACACCATAGAATGACTATGAGTCTCATAAGAATGCACTCAACTAGTCATCACTCCTGTGTTTTCATAAGAAAAAACAGTGTTAGAGTCCAAGAG",
                                "+",
                                "AAAAA6EEEEE/EEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
        }

        //void TearDown() override{
        //    delete rp.left;
        //    delete rp.right;
        //}

        fqlib::ReadPair rp;
};

TEST_F(ReadPairTest, merge){
    fqlib::Read* mr = rp.merge();
    EXPECT_EQ(mr->seq.seqStr, "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTTT");
    EXPECT_NE(mr->quality, "");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

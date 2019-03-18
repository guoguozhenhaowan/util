#include "seq.h"
#include "gtest/gtest.h"

// create SeqTest class
class SeqTest : public ::testing::Test{
    protected:
        void SetUp() override{
        }
        
        fqlib::Seq s1;
        fqlib::Seq s2 = {"ATCGGCAT"};
        fqlib::Seq s3 = {"GTACGTNU"};
};

// test default constructor
TEST_F(SeqTest, isEmptyInitially){
    EXPECT_EQ(s1.length(), 0);
}

// test reverseComplement function
TEST_F(SeqTest, reverseComplement){
    EXPECT_EQ(s2.reverseComplement().seqStr, "ATGCCGAT");
    EXPECT_EQ(s3.reverseComplement().seqStr, "NNACGTAC");
    fqlib::Seq s4 = ~s3;
    EXPECT_EQ(s4.seqStr, "NNACGTAC");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return ::RUN_ALL_TESTS();
}

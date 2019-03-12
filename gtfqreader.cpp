#include "fqreader.h"
#include "gtest/gtest.h"

class FastqReaderTest : public testing::Test{
    protected:
        void SetUp() override{
        }

        FqReader reader1 = {"./testdata/R1.fq"};
        FqReader reader2 = {"./testdata/R1.fq.gz"};
};

TEST_F(FastqReaderTest, read){
    EXPECT_FALSE(reader1.isZipped());
    EXPECT_TRUE(reader2.isZipped());
    size_t ttbt = 0, rdbt = 0;
    reader1.getBytes(rdbt, ttbt);
    EXPECT_EQ(rdbt, 3042);
    EXPECT_EQ(ttbt, 3042);
    reader2.getBytes(rdbt, ttbt);
    EXPECT_EQ(rdbt, 373);
    EXPECT_EQ(ttbt, 373);

    Read* r1 = NULL;
    Read* r2 = NULL;
    while(true){
        r1 = reader1.read();
        r2 = reader2.read();
        if(r1 && r2){
            EXPECT_EQ(r1->seq.seqStr, r2->seq.seqStr);
            //std::cout << "r1:" << r1->seq.seqStr << "\n";
            //std::cout << "r2:" << r2->seq.seqStr << "\n";
        }else{
            break;
        }
    }
    EXPECT_TRUE(reader1.eof());
    EXPECT_TRUE(reader2.eof());
    EXPECT_FALSE(reader1.hasNoLineBreakAtEnd());
    EXPECT_TRUE(reader2.hasNoLineBreakAtEnd());

    delete r1;
    delete r2;
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

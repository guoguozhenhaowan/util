#include "gtest/gtest.h"
#include "evaluator.h"

class EvaluatorTest : public testing::Test{
    void SetUp() override{
    }
    protected:
        Evaluator e = {"./testdata/R1.adaptor.fq.gz", 0};
};

TEST_F(EvaluatorTest, All){
    EXPECT_EQ(e.getAdapter(), "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    EXPECT_TRUE(e.isIlluminaAdapter());
    EXPECT_TRUE(e.isTwoColorSystem());
    EXPECT_EQ(e.getReadLen(), 150);
    EXPECT_TRUE(std::abs((long)e.getReadNum() - 300000) < 1000);
    std::cout << "estimated readnum: " << e.getReadNum() << std::endl;
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

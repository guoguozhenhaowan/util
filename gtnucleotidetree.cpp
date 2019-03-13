#include "nucleotidetree.h"
#include "gtest/gtest.h"

class NucleotideTreeTest : public testing::Test{
    protected:
        void SetUp() override{
            for(int i = 0; i < 100; ++i){
                ntree.addSeq("AAAATTTT");
                ntree.addSeq("AAAATTTTGGGG");
                ntree.addSeq("AAAATTTTGGGGCCCC");
                ntree.addSeq("AAAATTTTGGGGCCAA");
            }
            for(int i = 0; i < 49; ++i){
                ptree.addSeq("AAAATTTTGGGGCCCC");
            }
            ntree.addSeq("AAAATTTTGGGACCCC");
        }

        bool reachedLeaf = true;
        NucleotideTree ntree;
        NucleotideTree ptree;
};

TEST_F(NucleotideTreeTest, getDominantPath){
    std::string path = ntree.getDominantPath(reachedLeaf);
    EXPECT_EQ(path, "AAAATTTTGGGGCC");
    EXPECT_FALSE(reachedLeaf);
    path = ptree.getDominantPath(reachedLeaf);
    EXPECT_EQ(path, "");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

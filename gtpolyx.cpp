#include "polyx.h"
#include "gtest/gtest.h"

TEST(PolyX, test){
    fqlib::Read r("@name",
        "ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    fqlib::PolyX::trimPolyX(&r, NULL, 10);
    EXPECT_EQ(r.seq.seqStr, "ATTTT");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include "gtest/gtest.h"
#include "writer.h"
#include "fqreader.h"

TEST(Writer, write){
    std::string str = "fuck shit";
    char cstr[] = "hello kitty";
    std::string outf = "./plain.txt";
    std::string zipf = "./hello.gz";
    Writer fp = {outf};
    Writer fz = {zipf};
    fp.init();
    fp.writeLine(str);
    fp.write(cstr, 11);
    fp.writeString(str);
    fp.close();

    fz.init();
    fz.write(cstr, 10);
    fz.writeLine(str);
    fz.writeString(str);
    fz.close();

    FqReader fqr = {"./hello.gz"};
    std::string sget = fqr.getLine();
    EXPECT_EQ(sget, "hello kittfuck shit");

    FqReader fqr2 = {"./plain.txt"};
    sget = fqr2.getLine();
    EXPECT_EQ(sget, "fuck shit");
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

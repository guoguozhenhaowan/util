#include "util.h"
#include "gtest/gtest.h"

// test util::starts_with
TEST(util, starts_with){
    EXPECT_TRUE(util::starts_with("fuckshit", "fuck"));
    EXPECT_TRUE(util::starts_with(" hellokitty", " "));
    EXPECT_TRUE(util::starts_with(" ", ""));
    EXPECT_FALSE(util::starts_with("h", "hello"));
    EXPECT_FALSE(util::starts_with("hello", "ello"));
}

// test util::ends_with
TEST(util, ends_with){
    EXPECT_TRUE(util::ends_with("fuckshit", "it"));
    EXPECT_TRUE(util::ends_with("hellokit ", " "));
    EXPECT_TRUE(util::ends_with(" ", ""));
    EXPECT_FALSE(util::ends_with("h", "hello"));
    EXPECT_FALSE(util::ends_with("shit", "a"));
}

// test util::split
TEST(util, split){
    std::string s1("fuck shit    jjburst kao");
    std::string s2("\tkello\tshit\fburst\r");
    std::string s3(" k j\t k ");
    std::string s4(" ");
    std::string s5;
    std::vector<std::string> v1, v2, v3, v4, v5;
     
    util::split(s1, v1, " ");
    EXPECT_EQ(v1[0], "fuck");
    EXPECT_EQ(v1[1], "shit");
    EXPECT_EQ(v1[2], "jjburst");
    EXPECT_EQ(v1[3], "kao");
    
    util::split(s2, v2, "\t\f\r");
    EXPECT_EQ(v2[0], "kello");
    EXPECT_EQ(v2[1], "shit");
    EXPECT_EQ(v2[2], "burst");
    
    util::split(s3, v3, " \t");
    EXPECT_EQ(v3[0], "k");
    EXPECT_EQ(v3[1], "j");
    EXPECT_EQ(v3[0], v3[2]);

    util::split(s4, v4, " ");
    EXPECT_EQ(v4.size(), 0);

    util::split(s5, v5, " ");
    EXPECT_EQ(v5.size(), 0);
}

// test util::strip
TEST(util, strip){
    std::string s1("\nfuck shit    jjburst kao\t");
    std::string s2(" \t\nkello     shit burst");
    std::string s3(" \r\rk j k \f");
    std::string s4(" ");
    std::string s5;

    EXPECT_EQ(util::strip(s1), "fuck shit    jjburst kao");
    EXPECT_EQ(util::strip(s2), "kello     shit burst");
    EXPECT_EQ(util::strip(s3), "k j k");
    EXPECT_EQ(util::strip(s4), "");
    EXPECT_EQ(util::strip(s5), "");
}

// test util::rstrip
TEST(util, rstrip){
    std::string s1("fuck shit    jjburst kao ");
    std::string s2(" kello     shit burst\t \n ");
    std::string s3(" k j k\r");
    std::string s4(" ");
    std::string s5;

    EXPECT_EQ(util::rstrip(s1), "fuck shit    jjburst kao");
    EXPECT_EQ(util::rstrip(s2), " kello     shit burst");
    EXPECT_EQ(util::rstrip(s3), " k j k");
    EXPECT_EQ(util::rstrip(s4), "");
    EXPECT_EQ(util::rstrip(s5), s5);
}

// test util::lstrip
TEST(util, lstrip){
    std::string s1("fuck shit    jjburst kao");
    std::string s2(" kello     shit burst");
    std::string s3("\n\tk j k ");
    std::string s4(" ");
    std::string s5;

    EXPECT_EQ(util::lstrip(s1), s1);
    EXPECT_EQ(util::lstrip(s2), "kello     shit burst");
    EXPECT_EQ(util::lstrip(s3), "k j k ");
    EXPECT_EQ(util::lstrip(s4), "");
    EXPECT_EQ(util::lstrip(s5), s5);
}

// test util::replace
TEST(util, replace){
    std::string s1("fuck shit    jjburst kao");
    std::string s2(" kello    shit burst");
    std::string s3("\n\tk j k ");
    std::string s4(" ");
    std::string s5;

    EXPECT_EQ(util::replace(s1, "s", "x"), "fuck xhit    jjburxt kao");
    EXPECT_EQ(util::replace(s2, " ", "-"), "-kello----shit-burst");
    EXPECT_EQ(util::replace(s3, "k j", "fuckshit"), "\n\tfuckshit k ");
    EXPECT_EQ(util::replace(s4, " ", "k"), "k");
    EXPECT_EQ(util::replace(s5, " ", "u"), "");
}

// test util::basename
TEST(util, basename){
    std::string s1("fuck");
    std::string s2(".///fuck");
    std::string s3("./");
    std::string s4("/fuck//shit/hello///");
    std::string s5("/fuck/shit/kao");
    std::string s6("/fuck///");
    std::string s7(" ");

    EXPECT_EQ(util::basename(s1), s1);
    EXPECT_EQ(util::basename(s2), "fuck");
    EXPECT_EQ(util::basename(s3), ".");
    EXPECT_EQ(util::basename(s4), "hello");
    EXPECT_EQ(util::basename(s5), "kao");
    EXPECT_EQ(util::basename(s6), "fuck");
    EXPECT_EQ(util::basename(s7), "");
}

// test util::dirname
TEST(util, dirname){
    std::string s1("fuck");
    std::string s2(".///fuck");
    std::string s3("./");
    std::string s4("/fuck//shit/hello///");
    std::string s5("/fuck/shit/kao");
    std::string s6("/fuck///");
    std::string s7(" ");

    EXPECT_EQ(util::dirname(s1), "./");
    EXPECT_EQ(util::dirname(s2), ".///");
    EXPECT_EQ(util::dirname(s3), "./");
    EXPECT_EQ(util::dirname(s4), "/fuck//shit/");
    EXPECT_EQ(util::dirname(s5), "/fuck/shit/");
    EXPECT_EQ(util::dirname(s6), "/");
    EXPECT_EQ(util::dirname(s7), "./");
}

// test util::joinpath
TEST(util, joinpath){
    EXPECT_EQ(util::joinpath("/fuck", "hello"), "/fuck/hello");
    EXPECT_EQ(util::joinpath("./", "fuck"), ".//fuck");
}

// test util::isfile
TEST(util, isfile){
    EXPECT_TRUE(util::isfile("./util.h"));
    EXPECT_FALSE(util::isfile("./fuck"));
    EXPECT_FALSE(util::isfile("/Users/wood/"));
}

// test util::isdir
TEST(util, isdir){
    EXPECT_TRUE(util::isdir("./"));
    EXPECT_TRUE(util::isdir("/Users/wood/Documents/"));
    EXPECT_FALSE(util::isdir("./a.out"));
    EXPECT_FALSE(util::isdir("./fuck"));
}

// test util::makedir
TEST(util, makedir){
    EXPECT_TRUE(util::makedir("./"));
    EXPECT_TRUE(util::makedir("./fuck/shit/jj"));
    EXPECT_TRUE(util::makedir("~/hello/litty/jj"));
    EXPECT_FALSE(util::makedir("/fuck/shit"));
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include "duplicate.h"
#include "fqreader.h"
#include "gtest/gtest.h"

class DuplicateTest : public testing::Test{
    protected:
        void SetUp() override{
            opt.duplicate.enabled = true;
            opt.duplicate.keylen = 10;    //should be less than 16
            opt.duplicate.histSize = 32;
            hist = new int[opt.duplicate.histSize];
            memset(hist, 0, sizeof(int) * opt.duplicate.histSize);
            meanGC = new double[opt.duplicate.histSize];
            memset(meanGC, 0, sizeof(double) * opt.duplicate.histSize);
        }

        void TearDown() override{
            delete[] hist;
            delete[] meanGC;
        }

        fqlib::Options opt = {};
        fqlib::FqReader fqrse = {"./testdata/R1.adaptor.fq.gz"};
        fqlib::FqReaderPair fqrpe = {"./testdata/R1.adaptor.fq.gz", "./testdata/R2.adaptor.fq.gz"};
        int* hist;
        double* meanGC;
};

TEST_F(DuplicateTest, SE){
    fqlib::Duplicate dup(&opt);
    fqlib::Read* r;
    while((r = fqrse.read()) != NULL){
        dup.statRead(r);
    }
    double dpr = dup.statAll(hist, meanGC, opt.duplicate.histSize);
    for(int i = 1; i < opt.duplicate.histSize; ++i){
        std::cout << i << ":" << hist[i] << "\t" << meanGC[i] << "\n";
    }
    std::cout << "dup ratio: " << dpr << "\n";
}

TEST_F(DuplicateTest, PE){
    fqlib::Duplicate dup(&opt);
    fqlib::ReadPair* rp;
    while((rp = fqrpe.read()) != NULL){
        dup.statPair(rp->left, rp->right);
    }

    double dpr = dup.statAll(hist, meanGC, opt.duplicate.histSize);
    for(int i = 1; i < opt.duplicate.histSize; ++i){
        std::cout << i << ":" << hist[i] << "\t" << meanGC[i] << "\n";
    }
    std::cout << "dup ratio: " << dpr << "\n";
}

int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

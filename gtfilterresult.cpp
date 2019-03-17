#include "filterresult.h"
#include "gtest/gtest.h"

TEST(FilterResult, Test){
    FilterOpt *popt = new FilterOpt;
    popt->baseCorrection = true;
    popt->complexityFilter = true;
    popt->filterLongRead = true;
    popt->filterShortRead = true;
    popt->trimAdapter = true;
    popt->maxReadLen = 300;
    popt->minReadLen = 30;

    FilterResult ser = {popt, false};
    FilterResult per = {popt, true};

    ser.addAdapterTrimmed("ATCGATA");
    ser.addCorrection('c', 'a');
    ser.addCorrection('t', 'c');
    ser.addFilterResult(compar::FAIL_OVERLAP);
    ser.addFilterResult(compar::PASS_FILTER);
    ser.incCorrectedReads(10);
    std::ofstream fwse("./se.json");
    ser.reportAdaptersJsonSummary(fwse, "", "ATCCAT");
    ser.reportJsonBasic(fwse, "");
    fwse.close();
    fwse.open("./se.html");
    ser.reportAdaptersHtmlSummary(fwse, 30);
    ser.reportHtmlBasic(fwse, 30, 300);
    fwse.close();

    per.addAdapterTrimmed("ATCGATA", "CCCCCC");
    per.addCorrection('c', 'a');
    per.addCorrection('t', 'c');
    per.addFilterResult(compar::FAIL_OVERLAP);
    per.addFilterResult(compar::PASS_FILTER);
    per.incCorrectedReads(10);
    std::ofstream fwpe("./pe.json");
    per.reportAdaptersJsonSummary(fwpe, "", "ATCCAT", "TTCAT");
    per.reportJsonBasic(fwpe, "");
    fwpe.close();
    fwpe.open("./pe.html");
    per.reportAdaptersHtmlSummary(fwpe, 30);
    per.reportHtmlBasic(fwpe, 30, 300);
    fwpe.close();

}
    
int main(int argc, char** argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
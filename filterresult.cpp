#include "filterresult.h"

namespace fqlib{
    FilterResult::FilterResult(Options* opt, bool paired){
        mOptions = opt;
        mPaired = paired;
        mTrimmedAdapterBases = 0;
        mTrimmedAdapterReads = 0;
        for(int i = 0; i < compar::FILTER_RESULT_TYPES; ++i){
            mFilterReadStats[i] = 0;
        }
        mCorrectionMatrix = new size_t[64];
        std::memset(mCorrectionMatrix, 0, sizeof(size_t) * 64);
        mCorrectedReads = 0;
        mSummarized = false;
    }
    
    FilterResult::~FilterResult(){
        delete mCorrectionMatrix;
    }
    
    void FilterResult::addFilterResult(int result){
        if(result < compar::PASS_FILTER || result >= compar::FILTER_RESULT_TYPES){
            return;
        }
        if(mPaired){
            mFilterReadStats[result] += 2;
        }else{
            mFilterReadStats[result] += 1;
        }
    }
    
    FilterResult* FilterResult::merge(std::vector<FilterResult*>& list){
        if(list.size() == 0){
            return NULL;
        }
        FilterResult* result = new FilterResult(list[0]->mOptions, list[0]->mPaired);
        size_t* target = result->getFilterReadStats();
        size_t* correction = result->getCorrectionMatrix();
        for(int i = 0; i < list.size(); ++i){
            // update mFilterReadStats 
            size_t* curStats = list[i]->getFilterReadStats();
            for(int j = 0; j < compar::FILTER_RESULT_TYPES; ++j){
                target[j] += curStats[j];
            }
            // update mCorrectionMatrix
            size_t* curCorr = list[i]->getCorrectionMatrix();
            for(int p = 0; p < 64; ++p){
                correction[p] += curCorr[p];
            }
            // update mTrimmedAdapterReads/bases
            result->mTrimmedAdapterReads += list[i]->mTrimmedAdapterReads;
            result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;
    
            // update read counting
            result->mCorrectedReads += list[i]->mCorrectedReads;
    
            // merge adapter stats
            for(auto& e: list[i]->mAdapter1Count){
                if(result->mAdapter1Count.count(e.first) > 0){
                    ++result->mAdapter1Count[e.first];
                }else{
                    result->mAdapter1Count[e.first] = 0;
                }
            }
    
            for(auto& e: list[i]->mAdapter2Count){
                if(result->mAdapter2Count.count(e.first) > 0){
                    ++result->mAdapter2Count[e.first];
                }else{
                    result->mAdapter2Count[e.first] = 0;
                }
            }
        }
        return result;
    }
    
    void FilterResult::summary(bool force){
        if(mSummarized && !force){
            return;
        }
        mCorrectedBases = 0;
        for(int p = 0; p < 64; ++p){
            mCorrectedBases += mCorrectionMatrix[p];
        }
        mSummarized = true;
    }
    
    size_t FilterResult::getTotalCorrectedBases(){
        if(!mSummarized){
            summary();
        }
        return mCorrectedBases;
    }
    
    void FilterResult::addCorrection(char from, char to){
        int f = from & 0x07;
        int t = from & 0x07;
        ++mCorrectionMatrix[f * 8 +t];
    }
    
    void FilterResult::incCorrectedReads(int count){
        mCorrectedReads += count;
    }
    
    size_t FilterResult::getCorrectionNum(char from, char to){
        int f = from & 0x07;
        int t = to & 0x07;
        return mCorrectionMatrix[f * 8 + t];
    }
    
    void FilterResult::addAdapterTrimmed(const std::string& adapter, bool isR2){
        if(adapter.empty()){
            return;
        }
        ++mTrimmedAdapterReads;
        mTrimmedAdapterBases += adapter.length();
        if(!isR2){
            if(mAdapter1Count.count(adapter) > 0){
                ++mAdapter1Count[adapter];
            }else{
                mAdapter1Count[adapter] = 1;
            }
        }else{
            if(mAdapter2Count.count(adapter) > 0){
                ++mAdapter2Count[adapter];
            }else{
                mAdapter2Count[adapter] = 1;
            }
        }
    }
    
    void FilterResult::addAdapterTrimmed(const std::string& adapterR1, const std::string& adapterR2){
        mTrimmedAdapterReads += 2;
        mTrimmedAdapterBases += adapterR1.length() + adapterR2.length();
        if(!adapterR1.empty()){
            if(mAdapter1Count.count(adapterR1) > 0){
                ++mAdapter1Count[adapterR1];
            }else{
                mAdapter1Count[adapterR1] = 1;
            }
        }
    
        if(!adapterR2.empty()){
            if(mAdapter2Count.count(adapterR2) > 0){
                ++mAdapter2Count[adapterR2];
            }else{
                mAdapter2Count[adapterR2] = 1;
            }
        }
    }
    
    std::ostream& operator<<(std::ostream& os, FilterResult* re){
        Options* mOptions = re->mOptions;
        os << "reads passed filter: " << re->mFilterReadStats[compar::PASS_FILTER] << "\n";
        os << "reads failed due to low quality: " << re->mFilterReadStats[compar::FAIL_QUALITY] << "\n";
        os << "reads failed due to too many N: " << re->mFilterReadStats[compar::FAIL_N_BASE] << "\n";
        if(mOptions->lengthFilter.enabled){
            os << "reads failed due to too short: " << re->mFilterReadStats[compar::FAIL_LENGTH] << "\n";
            if(mOptions->lengthFilter.maxReadLength > 0){
                os << "reads failed due to too long: " << re->mFilterReadStats[compar::FAIL_TOO_LONG] << "\n";
            }
        }
        if(mOptions->complexityFilter.enabled){
            os << "reads failed due to low complexity: " << re->mFilterReadStats[compar::FAIL_COMPLEXITY] << "\n";
        }
        if(mOptions->adapter.enableTriming){
            os << "reads with adapter trimmed: " << re->mTrimmedAdapterReads << "\n";
            os << "bases trimmed due to adapters: " << re->mTrimmedAdapterBases << "\n";
        }
        if(mOptions->correction.enabled){
            os << "reads corrected by overlap analysis: " << re->mCorrectedReads << "\n";
            os << "bases corrected by overlap analysis: " << re->getTotalCorrectedBases() << "\n";
        }
        return os;
    }
    
    void FilterResult::reportJsonBasic(std::ofstream& ofs, const std::string& padding){
        ofs << "{\n";
        jsonutil::writeRecord(ofs, padding, "passed_filter_reads", mFilterReadStats[compar::PASS_FILTER]);
        jsonutil::writeRecord(ofs, padding, "low_quality_reads", mFilterReadStats[compar::FAIL_QUALITY]);
        jsonutil::writeRecord(ofs, padding, "too_many_N_reads", mFilterReadStats[compar::FAIL_N_BASE]);
        if(mOptions->correction.enabled){
            jsonutil::writeRecord(ofs, padding, "corrected_reads", mCorrectedReads);
            jsonutil::writeRecord(ofs, padding, "corrected_bases", getTotalCorrectedBases());
        }
        if(mOptions->complexityFilter.enabled){
            jsonutil::writeRecord(ofs, padding, "low_complexity_reads", mFilterReadStats[compar::FAIL_COMPLEXITY]);
        }
        if(mOptions->lengthFilter.enabled){
            jsonutil::writeRecord(ofs, padding, "too_short_reads", mFilterReadStats[compar::FAIL_LENGTH]);
            if(mOptions->lengthFilter.maxReadLength > 0){
                jsonutil::writeRecord(ofs, padding, "too_long_reads", mFilterReadStats[compar::FAIL_TOO_LONG]);
            }
        }
        ofs << padding << "}," << std::endl;
    }
    
    void FilterResult::reportHtmlBasic(std::ofstream& ofs, size_t totalReads, size_t totalBases){
        ofs << "<table class='summary_table'>\n";
        htmlutil::outputTableRow(ofs, "reads passed filters:", htmlutil::formatNumber(mFilterReadStats[compar::PASS_FILTER]) + " (" + std::to_string(mFilterReadStats[compar::PASS_FILTER] * 100.0 / totalBases) + "%)");
        htmlutil::outputTableRow(ofs, "low_quality_reads", htmlutil::formatNumber(mFilterReadStats[compar::FAIL_QUALITY]) + " (" + std::to_string(mFilterReadStats[compar::FAIL_QUALITY] * 100.0 / totalBases) + "%)");
        htmlutil::outputTableRow(ofs, "too_many_N_reads", htmlutil::formatNumber(mFilterReadStats[compar::FAIL_N_BASE]) + " (" + std::to_string(mFilterReadStats[compar::FAIL_N_BASE] * 100.0 / totalBases) + "%)");
        if(mOptions->correction.enabled){
            htmlutil::outputTableRow(ofs, "corrected_reads", htmlutil::formatNumber(mCorrectedReads) + " (" + std::to_string(mCorrectedReads * 100.0 / totalReads) + "%)");
            htmlutil::outputTableRow(ofs, "corrected_bases", htmlutil::formatNumber(getTotalCorrectedBases()) + " (" + std::to_string(getTotalCorrectedBases()) + " (" + std::to_string(getTotalCorrectedBases() * 100.0 / totalBases) + "%)");
        }
        if(mOptions->complexityFilter.enabled){
            htmlutil::outputTableRow(ofs, "low_complexity_reads", htmlutil::formatNumber(mFilterReadStats[compar::FAIL_COMPLEXITY]) + " (" + std::to_string(mFilterReadStats[compar::FAIL_COMPLEXITY] * 100.0 / totalReads) + "%)");
        }
        if(mOptions->lengthFilter.enabled){
            htmlutil::outputTableRow(ofs, "too_short_reads", htmlutil::formatNumber(mFilterReadStats[compar::FAIL_LENGTH]) + " (" + std::to_string(mFilterReadStats[compar::FAIL_LENGTH] * 100.0 / totalReads) + "%)");
            if(mOptions->lengthFilter.maxReadLength > 0){
                htmlutil::outputTableRow(ofs, "too_long_reads", htmlutil::formatNumber(mFilterReadStats[compar::FAIL_TOO_LONG]) + " (" + std::to_string(mFilterReadStats[compar::FAIL_TOO_LONG] * 100.0 /totalReads) + "%)");
            }
        }
    }
    
    void FilterResult::reportAdaptersJsonDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts){
        size_t totalAdapters = 0;
        for(auto& e: adapterCounts){
            totalAdapters += e.second;
        }
        if(totalAdapters == 0){
            return;
        }
        const double reportThreshold = 0.01;
        const double dTotalAdapters = (double)totalAdapters;
        bool firstItem = true;
        size_t reported = 0;
        for(auto& e: adapterCounts){
            if(e.second / dTotalAdapters < reportThreshold){
                continue;
            }
            if(!firstItem){
                ofs << ", ";
            }else{
                firstItem = false;
            }
            ofs << "\"" << e.first << "\":" << e.second;
            reported += e.second;
        }
        size_t unreported = totalAdapters - reported;
        if(unreported > 0){
            if(!firstItem){
                ofs << ", ";
            }else{
                ofs << "\"" << "others" << "\":" << unreported;
            }
        }
    }
    
    void FilterResult::reportAdaptersHtmlDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts, size_t totalBases){
        size_t totalAdapters = 0;
        size_t totalAdaptersBases = 0;
        for(auto& e: adapterCounts){
            totalAdapters += e.second;
            totalAdaptersBases += e.first.length();
        }
        double frac = (double)totalAdaptersBases / (double)totalBases;
        if(mPaired){
            frac *= 2.0;
        }
        if(frac < 0.01){
            ofs << "<div class='sub_section_tips'>The input has little adapter percentage (~" << std::to_string(frac * 100.0) << "%), probably it's trimmed before.</div>\n";
        }
        if(totalAdapters == 0){
            return;
        }
    
        ofs << "<table class='summary_table'>\n";
        ofs << "<tr><td class='adapter_col' style='font-size:14px;color:#ffffff;background:#556699'>" << "Sequence" << "</td><td class='col2' style='font-size:14px;color:#ffffff;background:#556699'>" << "Occurrences" << "</td    ></tr>\n";    
        const double reportThreshold = 0.01;
        const double dTotalAdapters = (double)totalAdapters;
        size_t reported = 0;
        for(auto& e: adapterCounts){
            if(e.second / dTotalAdapters < reportThreshold){
                continue;
            }
            htmlutil::outputTableRow(ofs, e.first, e.second);
            reported += e.second;
        }
        size_t unreported = totalAdapters - reported;
        if(unreported > 0){
            std::string tag = "other adapter sequences";
            if(reported == 0){
                tag = "all adapter sequences";
            }
            htmlutil::outputTableRow(ofs, tag, unreported);
        }
        ofs << "</table>\n";
    }
    
    void FilterResult::reportAdaptersJsonSummary(std::ofstream& ofs, const std::string& padding, const std::string& adapterSeq1, const std::string& adapterSeq2){
        ofs << "{\n";
        jsonutil::writeRecord(ofs, padding, "adapter_trimmed_reads", mTrimmedAdapterReads);
        jsonutil::writeRecord(ofs, padding, "adapter_trimmed_bases", mTrimmedAdapterBases);
        jsonutil::writeRecord(ofs, padding, "read1_adapter_sequence", adapterSeq1);
        if(mPaired){
            jsonutil::writeRecord(ofs, padding, "read2_adapter_sequence", adapterSeq2);
        }
        ofs << padding << "\t" << "\"read1_adapter_counts\": " << "{";
        reportAdaptersJsonDetails(ofs, mAdapter1Count);
        ofs << "}";
        if(mPaired){
            ofs << ",\n";
            ofs << padding << "\t" << "\"read2_adapter_counts\": " << "{";
            reportAdaptersJsonDetails(ofs, mAdapter2Count);
            ofs << "}\n";
        }
        ofs << padding << "}," << std::endl;
    }
    
    void FilterResult::reportAdaptersHtmlSummary(std::ofstream& ofs, size_t totalBases){
        ofs << "<div class='subsection_title' onclick=showOrHide('read1_adapters')>Adapter or bad ligation of read1</div>\n";
        ofs << "<div id='read1_adapters'>\n";
        reportAdaptersHtmlDetails(ofs, mAdapter1Count, totalBases);
        ofs << "</div>\n";
        if(mPaired){
            ofs << "<div class='subsection_title' onclick=showOrHide('read2_adapters')>Adapter or bad ligation of read2</div>\n";
            ofs << "<div id='read2_adapters'>\n";
            reportAdaptersHtmlDetails(ofs, mAdapter2Count, totalBases);
            ofs << "</div>\n";
        }
    }
}

#include "filterresult.h"

FilterResult::FilterResult(Options* opt, bool paired){
    this->opt = opt;
    this->paired = paired;
    this->trimmedAdapterBases = 0;
    this->trimmedAdapterReads = 0;
    for(int i = 0; i < compar::FILTER_RESULT_TYPES; ++i){
        this->filterReadStats[i] = 0;
    }
    this->correctionMatrix = new size_t[64];
    std::memset(this->correctionMatrix, 0, sizeof(size_t) * 64);
    this->correctedReads = 0;
    this->summarized = false;
}

FilterResult::~FilterResult(){
    delete this->correctionMatrix;
}

void FilterResult::addFilterResult(int result){
    if(result < compar::PASS_FILTER || result >= compar::FILTER_RESULT_TYPES){
        return;
    }
    if(this->paired){
        this->filterReadStats[result] += 2;
    }else{
        this->filterReadStats[result] += 1;
    }
}

FilterResult* FilterResult::merge(std::vector<FilterResult*>& list){
    if(list.size() == 0){
        return NULL;
    }
    FilterResult* result = new FilterResult(list[0]->opt, list[0]->paired);
    size_t* target = result->getFilterReadStats();
    size_t* correction = result->getCorrectionMatrix();
    for(int i = 0; i < list.size(); ++i){
        // update filterReadStats 
        size_t* curStats = list[i]->getFilterReadStats();
        for(int j = 0; j < compar::FILTER_RESULT_TYPES; ++j){
            target[j] += curStats[j];
        }
        // update correctionMatrix
        size_t* curCorr = list[i]->getCorrectionMatrix();
        for(int p = 0; p < 64; ++p){
            correction[p] += curCorr[p];
        }
        // update trimmedAdapterReads/bases
        result->trimmedAdapterReads += list[i]->trimmedAdapterReads;
        result->trimmedAdapterBases += list[i]->trimmedAdapterBases;

        // update read counting
        result->correctedReads += list[i]->correctedReads;

        // merge adapter stats
        for(auto& e: list[i]->adapter1){
            if(result->adapter1.count(e.first) > 0){
                ++result->adapter1[e.first];
            }else{
                result->adapter1[e.first] = 0;
            }
        }

        for(auto& e: list[i]->adapter2){
            if(result->adapter2.count(e.first) > 0){
                ++result->adapter2[e.first];
            }else{
                result->adapter2[e.first] = 0;
            }
        }
    }
    return result;
}

void FilterResult::summary(bool force){
    if(this->summarized && !force){
        return;
    }
    this->correctedBases = 0;
    for(int p = 0; p < 64; ++p){
        this->correctedBases += this->correctionMatrix[p];
    }
    this->summarized = true;
}

size_t FilterResult::getTotalCorrectedBases(){
    if(!this->summarized){
        this->summary();
    }
    return this->correctedBases;
}

void FilterResult::addCorrection(char from, char to){
    int f = from & 0x07;
    int t = from & 0x07;
    ++this->correctionMatrix[f * 8 +t];
}

void FilterResult::incCorrectedReads(int count){
    this->correctedReads += count;
}

size_t FilterResult::getCorrectionNum(char from, char to){
    int f = from & 0x07;
    int t = to & 0x07;
    return this->correctionMatrix[f * 8 + t];
}

void FilterResult::addAdapterTrimmed(const std::string& adapter, bool isR2){
    if(adapter.empty()){
        return;
    }
    ++this->trimmedAdapterReads;
    this->trimmedAdapterBases += adapter.length();
    if(!isR2){
        if(this->adapter1.count(adapter) > 0){
            ++this->adapter1[adapter];
        }else{
            this->adapter1[adapter] = 1;
        }
    }else{
        if(this->adapter2.count(adapter) > 0){
            ++this->adapter2[adapter];
        }else{
            this->adapter2[adapter] = 1;
        }
    }
}

void FilterResult::addAdapterTrimmed(const std::string& adapterR1, const std::string& adapterR2){
    this->trimmedAdapterReads += 2;
    this->trimmedAdapterBases += adapterR1.length() + adapterR2.length();
    if(!adapterR1.empty()){
        if(this->adapter1.count(adapterR1) > 0){
            ++this->adapter1[adapterR1];
        }else{
            this->adapter1[adapterR1] = 1;
        }
    }

    if(!adapterR2.empty()){
        if(this->adapter2.count(adapterR2) > 0){
            ++this->adapter2[adapterR2];
        }else{
            this->adapter2[adapterR2] = 1;
        }
    }
}

std::ostream& operator<<(std::ostream& os, FilterResult* re){
    Options* opt = re->opt;
    os << "reads passed filter: " << re->filterReadStats[compar::PASS_FILTER] << "\n";
    os << "reads failed due to low quality: " << re->filterReadStats[compar::FAIL_QUALITY] << "\n";
    os << "reads failed due to too many N: " << re->filterReadStats[compar::FAIL_N_BASE] << "\n";
    if(opt->lengthFilter.enabled){
        os << "reads failed due to too short: " << re->filterReadStats[compar::FAIL_LENGTH] << "\n";
        if(opt->lengthFilter.maxReadLength > 0){
            os << "reads failed due to too long: " << re->filterReadStats[compar::FAIL_TOO_LONG] << "\n";
        }
    }
    if(opt->complexityFilter.enabled){
        os << "reads failed due to low complexity: " << re->filterReadStats[compar::FAIL_COMPLEXITY] << "\n";
    }
    if(opt->adapter.enableTriming){
        os << "reads with adapter trimmed: " << re->trimmedAdapterReads << "\n";
        os << "bases trimmed due to adapters: " << re->trimmedAdapterBases << "\n";
    }
    if(opt->correction.enabled){
        os << "reads corrected by overlap analysis: " << re->correctedReads << "\n";
        os << "bases corrected by overlap analysis: " << re->getTotalCorrectedBases() << "\n";
    }
    return os;
}

void FilterResult::reportJsonBasic(std::ofstream& ofs, const std::string& padding){
    ofs << "{\n";
    jsonutil::writeRecord(ofs, padding, "passed_filter_reads", this->filterReadStats[compar::PASS_FILTER]);
    jsonutil::writeRecord(ofs, padding, "low_quality_reads", this->filterReadStats[compar::FAIL_QUALITY]);
    jsonutil::writeRecord(ofs, padding, "too_many_N_reads", this->filterReadStats[compar::FAIL_N_BASE]);
    if(opt->correction.enabled){
        jsonutil::writeRecord(ofs, padding, "corrected_reads", this->correctedReads);
        jsonutil::writeRecord(ofs, padding, "corrected_bases", this->getTotalCorrectedBases());
    }
    if(opt->complexityFilter.enabled){
        jsonutil::writeRecord(ofs, padding, "low_complexity_reads", this->filterReadStats[compar::FAIL_COMPLEXITY]);
    }
    if(opt->lengthFilter.enabled){
        jsonutil::writeRecord(ofs, padding, "too_short_reads", this->filterReadStats[compar::FAIL_LENGTH]);
        if(opt->lengthFilter.maxReadLength > 0){
            jsonutil::writeRecord(ofs, padding, "too_long_reads", this->filterReadStats[compar::FAIL_TOO_LONG]);
        }
    }
    ofs << padding << "}," << std::endl;
}

void FilterResult::reportHtmlBasic(std::ofstream& ofs, size_t totalReads, size_t totalBases){
    ofs << "<table class='summary_table'>\n";
    htmlutil::outputTableRow(ofs, "reads passed filters:", htmlutil::formatNumber(this->filterReadStats[compar::PASS_FILTER]) + " (" + std::to_string(this->filterReadStats[compar::PASS_FILTER] * 100.0 / totalBases) + "%)");
    htmlutil::outputTableRow(ofs, "low_quality_reads", htmlutil::formatNumber(this->filterReadStats[compar::FAIL_QUALITY]) + " (" + std::to_string(this->filterReadStats[compar::FAIL_QUALITY] * 100.0 / totalBases) + "%)");
    htmlutil::outputTableRow(ofs, "too_many_N_reads", htmlutil::formatNumber(this->filterReadStats[compar::FAIL_N_BASE]) + " (" + std::to_string(this->filterReadStats[compar::FAIL_N_BASE] * 100.0 / totalBases) + "%)");
    if(opt->correction.enabled){
        htmlutil::outputTableRow(ofs, "corrected_reads", htmlutil::formatNumber(this->correctedReads) + " (" + std::to_string(this->correctedReads * 100.0 / totalReads) + "%)");
        htmlutil::outputTableRow(ofs, "corrected_bases", htmlutil::formatNumber(this->getTotalCorrectedBases()) + " (" + std::to_string(this->getTotalCorrectedBases()) + " (" + std::to_string(this->getTotalCorrectedBases() * 100.0 / totalBases) + "%)");
    }
    if(opt->complexityFilter.enabled){
        htmlutil::outputTableRow(ofs, "low_complexity_reads", htmlutil::formatNumber(this->filterReadStats[compar::FAIL_COMPLEXITY]) + " (" + std::to_string(this->filterReadStats[compar::FAIL_COMPLEXITY] * 100.0 / totalReads) + "%)");
    }
    if(opt->lengthFilter.enabled){
        htmlutil::outputTableRow(ofs, "too_short_reads", htmlutil::formatNumber(this->filterReadStats[compar::FAIL_LENGTH]) + " (" + std::to_string(this->filterReadStats[compar::FAIL_LENGTH] * 100.0 / totalReads) + "%)");
        if(opt->lengthFilter.maxReadLength > 0){
            htmlutil::outputTableRow(ofs, "too_long_reads", htmlutil::formatNumber(this->filterReadStats[compar::FAIL_TOO_LONG]) + " (" + std::to_string(this->filterReadStats[compar::FAIL_TOO_LONG] * 100.0 /totalReads) + "%)");
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
    if(this->paired){
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
    jsonutil::writeRecord(ofs, padding, "adapter_trimmed_reads", this->trimmedAdapterReads);
    jsonutil::writeRecord(ofs, padding, "adapter_trimmed_bases", this->trimmedAdapterBases);
    jsonutil::writeRecord(ofs, padding, "read1_adapter_sequence", adapterSeq1);
    if(this->paired){
        jsonutil::writeRecord(ofs, padding, "read2_adapter_sequence", adapterSeq2);
    }
    ofs << padding << "\t" << "\"read1_adapter_counts\": " << "{";
    this->reportAdaptersJsonDetails(ofs, this->adapter1);
    ofs << "}";
    if(this->paired){
        ofs << ",\n";
        ofs << padding << "\t" << "\"read2_adapter_counts\": " << "{";
        this->reportAdaptersJsonDetails(ofs, this->adapter2);
        ofs << "}\n";
    }
    ofs << padding << "}," << std::endl;
}

void FilterResult::reportAdaptersHtmlSummary(std::ofstream& ofs, size_t totalBases){
    ofs << "<div class='subsection_title' onclick=showOrHide('read1_adapters')>Adapter or bad ligation of read1</div>\n";
    ofs << "<div id='read1_adapters'>\n";
    this->reportAdaptersHtmlDetails(ofs, this->adapter1, totalBases);
    ofs << "</div>\n";
    if(this->paired){
        ofs << "<div class='subsection_title' onclick=showOrHide('read2_adapters')>Adapter or bad ligation of read2</div>\n";
        ofs << "<div id='read2_adapters'>\n";
        this->reportAdaptersHtmlDetails(ofs, this->adapter2, totalBases);
        ofs << "</div>\n";
    }
}

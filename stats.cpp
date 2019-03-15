#include "stats.h"

Stats::Stats(const int& estReadLen){
    this->reads = 0;
    this->bases = 0;
    this->evaluatedSeqLen = estReadLen;
    this->cycles = estReadLen;
    this->bufLen = estReadLen + 1024;
    this->q20Total = 0;
    this->q30Total = 0;
    this->summarized = false;
    this->kmerMax = 0;
    this->kmerMin = 0;
    this->kmerLen = 5;
    this->kmerBufLen = 2 << (this->kmerLen * 2);
    this->lengthSum = 0;
    this->overRepSampleFreq = 100;
}

void Stats::allocateRes(){
    for(int i = 0; i < 8; ++i){
        this->q20Bases[i] = 0;
        this->q30Bases[i] = 0;
        this->baseContents[i] = 0;
        this->cycleQ20Bases[i] = new size_t[this->bufLen];
        std::memset(this->cycleQ20Bases[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleQ30Bases[i] = new size_t[this->bufLen];
        std::memset(this->cycleQ30Bases[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleBaseContents[i] = new size_t[this->bufLen];
        std::memset(this->cycleBaseContents[i], 0, sizeof(size_t) * this->bufLen);
        this->cycleBaseQual[i] = new size_t[this->bufLen];
        std::memset(this->cycleBaseQual[i], 0, sizeof(size_t) * this->bufLen);
    }

    this->cycleTotalBase = new size_t[this->bufLen];
    std::memset(this->cycleTotalBase, 0, sizeof(size_t) * this->bufLen);
    this->cycleTotalQual = new size_t[this->bufLen];
    std::memset(this->cycleTotalQual, 0, sizeof(size_t) * this->bufLen);
    
    if(this->kmerLen){
        this->kmer = new size_t[this->kmerBufLen];
        std::memset(this->kmer, 0, sizeof(size_t) * this->kmerBufLen);
    }
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= this->bufLen){
        return;
    }
    size_t* newBuf = NULL;

    for(int i = 0; i < 8; ++i){
        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleQ20Bases[i], sizeof(size_t) * this->bufLen);
        delete this->cycleQ20Bases[i];
        this->cycleQ20Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleQ30Bases[i], sizeof(size_t) * this->bufLen);
        delete this->cycleQ30Bases[i];
        this->cycleQ30Bases[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleBaseContents[i], sizeof(size_t) * this->bufLen);
        delete this->cycleBaseContents[i];
        this->cycleBaseContents[i] = newBuf;

        newBuf = new size_t[newBufLen];
        std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
        std::memcpy(newBuf, this->cycleBaseQual[i], sizeof(size_t) * this->bufLen);
        delete this->cycleBaseQual[i];
        this->cycleBaseQual[i] = newBuf;
    }

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, this->cycleTotalBase, sizeof(size_t) * this->bufLen);
    delete this->cycleTotalBase;
    this->cycleTotalBase = newBuf;

    newBuf = new size_t[newBufLen];
    std::memset(newBuf, 0, sizeof(size_t) * newBufLen);
    memcpy(newBuf, this->cycleTotalQual, sizeof(size_t) * this->bufLen);
    delete this->cycleTotalQual;
    this->cycleTotalQual = newBuf;

    this->bufLen = newBufLen;
}

Stats::~Stats(){
    for(int i = 0; i < 8; ++i){
        delete cycleQ20Bases[i];
        cycleQ20Bases[i] = NULL;
        
        delete cycleQ30Bases[i];
        cycleQ30Bases[i] = NULL;
        
        delete cycleBaseContents[i];
        cycleBaseContents[i] = NULL;
        
        delete cycleBaseQual[i];
        cycleBaseQual[i] = NULL;
    }
    
    delete cycleTotalBase;
    delete cycleTotalQual;
    
    for(auto& e : this->qualityCurves){
        delete e.second;
    }
    
    for(auto& e : this->contentCurves){
        delete e.second;
    }
    
    if(this->kmerLen){
        delete this->kmer;
    }

    this->deleteOverRepSeqDist();
}

void Stats::summarize(bool forced){
    if(this->summarized && !forced){
        return;
    }

    // cycle and total bases
    int c = 0;
    for(c = 0; c < this->bufLen; ++c){
        this->bases += this->cycleTotalBase[c];
        if(this->cycleTotalBase[c] == 0){
            break;
        }
    }
    this->cycles = c;

    // Q20, Q30, base content
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < this->cycles; ++j){
            this->q20Bases[i] += this->cycleQ20Bases[i][j];
            this->q30Bases[i] += this->cycleQ30Bases[i][j];
            this->baseContents[i] += this->cycleBaseContents[i][j];
        }
        this->q20Total += this->q20Bases[i];
        this->q30Total += this->q30Bases[i];
    }

    // quality curve for mean qual
    double* meanQualCurve = new double[this->cycles];
    std::memset(meanQualCurve, 0, sizeof(double) * this->cycles);
    for(int i = 0; i < this->cycles; ++i){
        meanQualCurve[i] = (double)this->cycleTotalQual[i] / (double)this->cycleTotalBase[i];
    }
    this->qualityCurves["mean"] = meanQualCurve;

    // quality curves and base contents curves for different nucleotides
    char nucleotides[5] = {'A', 'T', 'C', 'G', 'N'};
    for(int i = 0; i < 5; ++i){
        char b = nucleotides[i] & 0x07;
        double* qualCurve = new double[this->cycles];
        std::memset(qualCurve, 0, sizeof(double) * this->cycles);
        double* contentCurve = new double[this->cycles];
        std::memset(contentCurve, 0, sizeof(double) * this->cycles);
        for(int j = 0; j < this->cycles; ++j){
            qualCurve[j] = (double)this->cycleBaseQual[b][j] / (double)this->cycleBaseContents[b][j];
            contentCurve[j] = (double)this->cycleBaseContents[b][j] / (double)this->cycleBaseContents[b][j];
        }
        this->qualityCurves[std::string(1, nucleotides[i])] = qualCurve;
        this->contentCurves[std::string(1, nucleotides[i])] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[this->cycles];
    std::memset(gcContentCurve, 0, sizeof(double) * this->cycles);
    char gIndex = 'G' & 0x07;
    char cIndex = 'C' & 0x07;
    for(int i = 0; i < this->cycles; ++i){
        gcContentCurve[i] = (double)(this->cycleBaseContents[gIndex][i] + this->cycleBaseContents[cIndex][i])/ (double)this->cycleTotalBase[i];
    }
    this->contentCurves["GC"] = gcContentCurve;

    if(this->kmerLen){
        // K-mer statistics
        this->kmerMin = this->kmer[0];
        this->kmerMax = this->kmer[0];
        for(int i = 1; i < this->kmerBufLen; ++i){
            this->kmerMin = std::min(this->kmer[i], this->kmerMin);
            this->kmerMax = std::max(this->kmer[i], this->kmerMax);
        }
    }

    this->summarized = true;
}

int Stats::getMeanLength(){
    if(this->reads == 0){
        return 0;
    }
    return this->lengthSum/this->reads;
}

void Stats::statRead(Read* r){
    int len = r->length();
    this->lengthSum += len;
    if(this->bufLen < len){
        this->extendBuffer(std::max(len + 100, (int)(len * 1.5)));
    }
    const std::string seqStr = r->seq.seqStr.c_str();
    const std::string qualStr = r->quality.c_str();
    int kmerVal = -1;
    for(int i = 0; i < len; ++i){
        char base = seqStr[i];
        char qual = qualStr[i];
        char bIndex = base & 0x07;
        const char q20 = '5';
        const char q30 = '?';
        if(qual > q20){
            ++(this->cycleQ20Bases[bIndex][i]);
        }
        if(qual > q30){
            ++this->cycleQ20Bases[bIndex][i];
            ++this->cycleQ30Bases[bIndex][i];
        }
        ++this->cycleBaseContents[bIndex][i];
        this->cycleBaseQual[bIndex][i] += (qual - 33);
        ++this->cycleTotalBase[i];
        this->cycleTotalQual[i] += (qual - 33);
        
        if(this->kmerLen){
            if(i < this->kmerLen - 1){
                continue;
            }
            kmerVal = Evaluator::seq2int(seqStr, i - this->kmerLen + 1, this->kmerLen, kmerVal);
            if(kmerVal >= 0){
                ++this->kmer[kmerVal];
            }
        }

        if(this->overRepSampleFreq){
            if(this->reads % this->overRepSampleFreq == 0){
                std::set<int> steps = {10, 20, 40, 100, std::min(150, this->evaluatedSeqLen - 2)};
                for(auto step: steps){
                    for(int i = 0; i < len - step; i += step){
                        std::string seq = seqStr.substr(i, step);
                        if(this->overRepSeq.count(seq) > 0){
                            ++this->overRepSeq[seq];
                        }else{
                            for(int p = i; p < seq.length() + i && p < this->evaluatedSeqLen; ++p){
                                if(this->overRepSeqDist.find(seq) == this->overRepSeqDist.end()){
                                    this->overRepSeqDist[seq] = new size_t[this->evaluatedSeqLen];
                                    std::memset(this->overRepSeqDist[seq], 0, sizeof(size_t) * this->evaluatedSeqLen);
                                }
                                ++this->overRepSeqDist[seq][p];
                            }
                        }
                    }
                }
            }
        }
    }
    ++this->reads;
}

int Stats::getCycles(){
    if(!this->summarized){
        this->summarize();
    }
    return this->cycles;
}

size_t Stats::getReads(){
    if(!this->summarized){
        this->summarize();
    }
    return this->reads;
}

size_t Stats::getBases(){
    if(!this->summarized){
        this->summarize();
    }
    return this->bases;
}

size_t Stats::getQ20(){
    if(!this->summarized){
        this->summarize();
    }
    return this->q20Total;
}

size_t Stats::getQ30(){
    if(!this->summarized){
        this->summarize();
    }
    return this->q30Total;
}

size_t Stats::getGCNumber(){
    if(!this->summarized){
        this->summarize();
    }
    return this->baseContents['G' & 0x07] + this->baseContents['C' & 0x07];
}

std::ostream& operator<<(std::ostream& os, const Stats& s){
    os << "total reads: " << s.reads << "\n";
    os << "total bases: " << s.bases << "\n";
    os << "Q20   bases: " << s.q20Total << "(" << (s.q20Total*100.0)/s.bases << "%)\n";
    os << "Q30   bases: " << s.q30Total << "(" << (s.q30Total*100.0)/s.bases << "%)\n";
    return os;
}

void Stats::reportJson(std::ofstream& ofs, std::string padding){
    // braces
    ofs << "{\n";
    // basic qc info 
    ofs << padding << "\t\"total_reads\": " << this->reads << ",\n";
    ofs << padding << "\t\"total_bases\": " << this->bases << ",\n";
    ofs << padding << "\t\"q20_bases\": " << this->q20Total << ",\n";
    ofs << padding << "\t\"q30_bases\": " << this->q30Total << ",\n";
    ofs << padding << "\t\"total_cycles\": " << this->cycles << ",\n";
    // quality curves
    std::string qualNames[5] = {"A", "T", "C", "G", "mean"};
    ofs << padding << "\t\"quality_curves\": {\n";
    for(int i = 0; i < 5; ++i){
        std::string name = qualNames[i];
        double* curve = this->qualityCurves[name];
        ofs << padding << "\t\t\"" << name << "\":[";
        for(int c = 0; c < this->cycles; ++c){
            ofs << curve[c];
            if(c != this->cycles - 1){
                ofs << ",";
            }
        }
        ofs << "]";
        if(i != 4){
            ofs << ",";
        }
        ofs << "\n";
    }
    ofs << padding << "\t},\n";
    //content curves
    std::string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    ofs << padding << "\t\"content_curves\":{\n";
    for(int i = 0; i < 6; ++i){
        std::string name = contentNames[i];
        double* curve = this->contentCurves[name];
        ofs << padding << "\t\t" << name << "\":[";
        for(int c = 0; c < this->cycles; ++c){
            ofs << curve[c];
            if(c != this->cycles - 1){
                ofs << ",";
            }
        }
        ofs << "]";
        if(i != 5){
            ofs << ",";
        }
        ofs << "\n";
    }
    ofs << padding << "\t},\n";
    //KMER counting(optional)
    if(this->kmerLen){
        ofs << padding << "\t" << "\"kmer_count\": {\n";
        for(size_t i = 0; i < (2 << (this->kmerLen * 2)); ++i){
            std::string seq = Evaluator::int2seq(i, this->kmerLen);
            ofs << padding << "\t\t\"" << seq << "\":" << this->kmer[i];
            if(i != (2 << this->kmerLen) - 1){
                ofs << ",";
            }
            if(i != (2 << (this->kmerLen * 2)) - 1){
                ofs << ",\n";
            }else{
                ofs << "\n";
            }
        }
        ofs << padding << "\t},\n";
    }
    // over represented seqs(optional)
    if(this->overRepSampleFreq){
        ofs << padding << "\t\"overrepresented_sequences\": {\n";
        bool first = true;
        for(auto& e : this->overRepSeq){
            if(!this->overRepPassed(e.first, e.second)){
                continue;
            }
            if(!first){
                ofs << ",\n";
            }else{
                first = false;
            }
            ofs << padding << "\t\t\"" << e.first << "\":" << e.second;
        }
        ofs << padding << "\t}\n";
    }
    ofs << padding << "},\n";
}

template<typename T>
std::string Stats::list2string(T* list, int size){
    std::stringstream ss;
    for(int i = 0; i < size; ++i){
        ss << list[i];
        if(i < size - 1){
            ss << ",";
        }
    }
    return ss.str();
}

template<typename T>
std::string Stats::list2string(T* list, int size, size_t* coords){
    std::stringstream ss;
    long start = 0;
    long end = 0;
    T total = 0;
    for(int i = 0; i < size; ++i){
        if(i > 0){
            start = coords[i - 1];
        }
        end = coords[i];
        total = 0;
        for(int j = start; j < end; ++j){
            total += list[j];
        }
        if(end == start){
            ss << "0";
        }else{
            ss << total /(end - start);
        }
        if(i < size - 1){
            ss << ",";
        }
    }
    return ss.str();
}

bool Stats::overRepPassed(const std::string& seq, size_t count){
    int s = this->overRepSampleFreq;
    switch(seq.length()){
        case 10:
            return s * count > 500;
        case 20:
            return s * count > 200;
        case 40:
            return s * count > 100;
        case 100:
            return s * count > 50;
        default:
            return s * count > 20;
    }
}

bool Stats::isLongRead(){
    return this->cycles > 300;
}

void Stats::reportHtml(std::ofstream& ofs, std::string filteringType, std::string readName){
    this->reportHtmlQuality(ofs, filteringType, readName);
    this->reportHtmlContents(ofs, filteringType, readName);
    if(this->kmerLen){
        this->reportHtmlKmer(ofs, filteringType, readName);
    }
    if(this->overRepSampleFreq){
        this->reportHtmlORA(ofs, filteringType, readName);
    }
}

void Stats::reportHtmlORA(std::ofstream& ofs, std::string filteringType, std::string readName){
    double dBases = this->bases;
    int displayed = 0;
    std::string subSection = filteringType + ": " + readName + ": overrepresented sequences";
    std::string divName = util::replace(subSection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";
    
    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('";
    ofs << divName << "')>" << subSection << "</a></div>\n";
    ofs << "<div id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Sampling rate: 1 / " << this->overRepSampleFreq << "</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr style='font-weight:bold;'><td>overrepresented sequence</td><td>count (% of bases)</td><td>distribution: cycle 1 ~ cycle " << this->evaluatedSeqLen << "</td></tr>"<<std::endl;
    int found = 0;
    for(auto& e: this->overRepSeq){
        std::string seq = e.first;
        size_t count = e.second;
        if(!this->overRepPassed(seq, count))
            continue;
        found++;
        double percent = (100.0 * count * seq.length() * this->overRepSampleFreq)/dBases;
        ofs << "<tr>";
        ofs << "<td width='400' style='word-break:break-all;font-size:8px;'>" << seq << "</td>";
        ofs << "<td width='200'>" << count << " (" << std::to_string(percent) <<"%)</td>";
        ofs << "<td width='250'><canvas id='" << divName << "_" << seq << "' width='240' height='20'></td>";
        ofs << "</tr>" << std::endl;
    }
    if(found == 0)
        ofs << "<tr><td style='text-align:center' colspan='3'>not found</td></tr>" << std::endl;
    ofs << "</table>\n";
    ofs << "</div>\n";

    // output the JS
    ofs << "<script language='javascript'>" << std::endl;
    ofs << "var seqlen = " << this->evaluatedSeqLen << ";" << std::endl;
    ofs << "var orp_dist = {" << std::endl;
    bool first = true;
    for(auto& e: this->overRepSeq){
        std::string seq = e.first;
        size_t count = e.second;
        if(!overRepPassed(seq, count))
            continue;
        if(!first) {
            ofs << "," << std::endl;
        } else
            first = false;
        ofs << "\t\"" << divName << "_" << seq << "\":[";
        for(int i=0; i< this->evaluatedSeqLen; ++i){
            if(i !=0 )
                ofs << ",";
            ofs << this->overRepSeqDist[seq][i];
        }
        ofs << "]";
    }
        ofs << "\n};" << std::endl;

    ofs << "for (seq in orp_dist) {"<< std::endl;
    ofs << "    var cvs = document.getElementById(seq);"<< std::endl;
    ofs << "    var ctx = cvs.getContext('2d'); "<< std::endl;
    ofs << "    var data = orp_dist[seq];"<< std::endl;
    ofs << "    var w = 240;"<< std::endl;
    ofs << "    var h = 20;"<< std::endl;
    ofs << "    ctx.fillStyle='#cccccc';"<< std::endl;
    ofs << "    ctx.fillRect(0, 0, w, h);"<< std::endl;
    ofs << "    ctx.fillStyle='#0000FF';"<< std::endl;
    ofs << "    var maxVal = 0;"<< std::endl;
    ofs << "    for(d=0; d<seqlen; d++) {"<< std::endl;
    ofs << "        if(data[d]>maxVal) maxVal = data[d];"<< std::endl;
    ofs << "    }"<< std::endl;
    ofs << "    var step = (seqlen-1) /  (w-1);"<< std::endl;
    ofs << "    for(x=0; x<w; x++){"<< std::endl;
    ofs << "        var target = step * x;"<< std::endl;
    ofs << "        var val = data[Math.floor(target)];"<< std::endl;
    ofs << "        var y = Math.floor((val / maxVal) * h);"<< std::endl;
    ofs << "        ctx.fillRect(x,h-1, 1, -y);"<< std::endl;
    ofs << "    }"<< std::endl;
    ofs << "}"<< std::endl;
    ofs << "</script>"<< std::endl;
}

void Stats::reportHtmlKmer(std::ofstream& ofs, std::string filteringType, std::string readName) {
    // KMER
    std::string subsection = filteringType + ": " + readName + ": KMER counting";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    ofs << "<table class='kmer_table' style='width:680px;'>\n";
    ofs << "<tr>";
    ofs << "<td></td>";
    // the heading row
    for(int h = 0; h < (2 << this->kmerLen); ++h){
        ofs << "<td style='color:#333333'>" << std::to_string(h + 1) << "</td>";
    }
    ofs << "</tr>\n";
    // content
    for(int i = 0; i < (2 << this->kmerLen); ++i){
        ofs << "<tr>";
        ofs << "<td style='color:#333333'>" << std::to_string(i + 1) << "</td>";
        for(int j = 0; j < (2 << this->kmerLen); ++j){
            ofs << this->makeKmerTD(i, j);
        }
        ofs << "</tr>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}

std::string Stats::makeKmerTD(int i, int j){
    std::string seq = Evaluator::int2seq(i * j, this->kmerLen);
    double meanBases = (double)(this->bases + 1) / this->kmerBufLen;
    double prop = this->kmer[i * j] / meanBases;
    double frac = 0.5;
    if(prop > 2.0){
        frac = (prop-2.0)/20.0 + 0.5;
    }
    else if(prop< 0.5){
        frac = prop;
    }

    frac = std::max(0.01, std::min(1.0, frac));
    int r = (1.0-frac) * 255;
    int g = r;
    int b = r;
    std::stringstream ss;
    ss << "<td style='background:#";
    if(r<16){
        ss << "0";
    }
    ss << std::hex <<r;
    if(g<16){
        ss << "0";
    }
    ss << std::hex <<g;
    if(b<16){
        ss << "0";
    }
    ss << std::hex << b;
    ss << std::dec << "' title='"<< seq << ": " << this->kmer[i * j] << "\n" << prop << " times as mean value'>";
    ss << seq << "</td>";
    return ss.str();
}

void Stats::reportHtmlQuality(std::ofstream& ofs, std::string filteringType, std::string readName) {
    // quality
    std::string subsection = filteringType + ": " + readName + ": quality";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[5] = {"A", "T", "C", "G", "mean"};
    std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    size_t *x = new size_t[this->cycles];
    int total = 0;
    if(!this->isLongRead()) {
        for(int i = 0; i < this->cycles; ++i){
            x[total] = i + 1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i < fullSampling && i < this->cycles; i++){
            x[total] = i + 1;
            total++;
        }
        // down sampling if it's too long
        if(this->cycles > fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= this->cycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != this->cycles){
                x[total] = this->cycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<5; b++) {
        std::string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + this->list2string(x, total) + "],";
        json_str += "y:[" + this->list2string(this->qualityCurves[base], total, x) + "],";
        json_str += "name: '" + base + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(this->isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'quality'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
}

void Stats::reportHtmlContents(std::ofstream& ofs, std::string filteringType, std::string readName) {

    // content
    std::string subsection = filteringType + ": " + readName + ": base contents";
    std::string divName = util::replace(subsection, " ", "_");
    divName = util::replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    size_t *x = new size_t[this->cycles];
    int total = 0;
    if(!this->isLongRead()) {
        for(int i=0; i<this->cycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<this->cycles; i++){
            x[total] = i+1;
             total++;
        }
        // down sampling if it's too long
        if(this->cycles > fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= this->cycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != this->cycles){
                x[total] = this->cycles;
                total++; 
             }
        } 
    }
    // four bases
    for (int b = 0; b<6; b++) {
        std::string base = alphabets[b];
        long count = 0;
        if(base.size()==1) {
            char b = base[0] & 0x07;
            count = this->baseContents[b];
        } else {
            count = this->baseContents['G' & 0x07] + this->baseContents['C' & 0x07] ;
        }
        std::string percentage = std::to_string((double)count * 100.0 / this->bases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        std::string name = base + "(" + percentage + "%)";

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(this->contentCurves[base], total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(this->isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
}

Stats* Stats::merge(std::vector<Stats*>& list, const int& estReadLen){
    if(list.size() == 0){
        return NULL;
    }
    
    int c = 0;
    int kl = 0;
    int smp = 0;
    for(int i = 0; i < list.size(); ++i){
        list[i]->summarize();
        c = std::max(c, list[i]->getCycles());
        kl = std::max(kl, list[i]->kmerLen);
        smp = std::max(smp, list[i]->overRepSampleFreq);
    }

    Stats* s = new Stats(estReadLen);
    s->setKmerLen(kl);
    s->setOverRepSampleFreq(smp);
    s->allocateRes();

    for(int i = 0; i < list.size(); ++i){
        int curCycles = list[i]->getCycles();
        s->reads += list[i]->reads;
        s->lengthSum += list[i]->lengthSum;
        for(int j = 0; j < 8; ++j){
            for(int k = 0; k < c && k < curCycles; ++k){
                s->cycleQ30Bases[i][j] += list[i]->cycleQ30Bases[i][j];
                s->cycleQ20Bases[i][j] += list[i]->cycleQ20Bases[i][j];
                s->cycleBaseContents[i][j] += list[i]->cycleBaseContents[i][j];
                s->cycleBaseQual[i][j] += list[i]->cycleBaseQual[i][j];
            }
        }

        for(int j = 0; j < c && j < curCycles; ++j){
            s->cycleTotalBase[j] += list[i]->cycleTotalBase[j];
            s->cycleTotalQual[j] += list[i]->cycleTotalQual[j];
        }

        if(s->kmerLen){
            for(int j = 0; j < s->kmerBufLen; ++j){
                s->kmer[i] += list[i]->kmer[i];
            }
        }

        if(s->overRepSampleFreq){
            for(auto& e: s->overRepSeq){
                s->overRepSeq[e.first] += list[i]->overRepSeq[e.first];
                for(int j = 0; j < s->evaluatedSeqLen; ++j){
                    if(s->overRepSeqDist.find(e.first) == s->overRepSeqDist.end()){
                        s->overRepSeqDist[e.first] = new size_t[s->evaluatedSeqLen];
                        std::memset(s->overRepSeqDist[e.first], 0, sizeof(size_t) * s->evaluatedSeqLen);
                    }
                    s->overRepSeqDist[e.first][j] += list[i]->overRepSeqDist[e.first][j];
                }
            }
        }
    }

    s->summarize();
    return s;
}

void Stats::deleteOverRepSeqDist(){
    for(auto& e: this->overRepSeq){
        delete this->overRepSeqDist[e.first];
        this->overRepSeqDist[e.first] = NULL;
    }
}

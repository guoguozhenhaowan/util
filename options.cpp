#include "options.h"

Options::Options(){
    this->in1 = "";
    this->in2 = "";
    this->out1 = "";
    this->out2 = "";
    this->reportTitle = "Fastq Report";
    this->thread = 1;
    this->compression = 2;
    this->phred64 = false;
    this->donotOverwrite = false;
    this->inputFromSTDIN = false;
    this->outputToSTDOUT = false;
    this->readsToProcess = 0;
    this->interleavedInput = false;
    this->insertSizeMax = 512;
    this->overlapDiffLimit = 5;
    this->verbose = false;
}

void Options::init(){
}

bool Options::isPaired(){
    return this->in2.length() > 0 || this->interleavedInput;
}

bool Options::adapterCutEnabled(){
    if(adapter.enableTriming){
        if(this->isPaired() || !adapter.inputAdapterSeqR1.empty()){
            return true;
        }
    }
    return false;
}

bool Options::validate(){
    // check in1/2 
    if(this->in1.empty()){
        if(!this->in2.empty()){
            util::error_exit("read2 input is specified by <in2>, but read1 input is not not specified by <in1>");
        }
        if(this->inputFromSTDIN){
            in1 = "/dev/stdin";
        }else{
            util::error_exit("read1 input should be specified by --in1, or enable --stdin if you want to read STDIN");
        }
    }else{
        util::isfile(this->in1);
    }
    if(!in2.empty()){
        util::isfile(this->in2);
    }
    
    // if output to STDOUT
    if(this->outputToSTDOUT){
        std::cerr << "Streaming uncompressed output to STDOUT...\n";
        if(!this->in1.empty() && !this->in2.empty()){
            std::cerr << "Enable interleaved output mode for paired-end input.\n";
        }
        if(!this->out1.empty()){
            std::cerr << "Ignore argument --out1 = " << this->out1 << "\n";
            this->out1 = "";
        }
        if(!this->out2.empty()){
            std::cerr << "Ignore argument --out2 = " << this->out2 << "\n";
            this->out2 = "";
        }
        if(this->split.enabled){
            std::cerr << "Ignore split mode\n";
            this->split.enabled = false;
        }
        std::cerr << std::endl;
    }
   
    // check in1/2
    if(this->in2.empty() && !this->interleavedInput && !this->out2.empty()){
        util::error_exit("read2 output is specified (--out2), but neighter read2 input is not specified (--in2), nor read1 is interleaved.");
    }

    if(!this->in2.empty() || this->interleavedInput){
        if(!this->out1.empty() && this->out2.empty()){
            util::error_exit("paired-end input, read1 output should be specified together with read2 output (--out2 needed) ");
        }
        if(this->out1.empty() && !this->out2.empty()){
            util::error_exit("paired-end input, read1 output should be specified (--out1 needed) together with read2 output ");
        }
    }
    if(!this->in2.empty() && this->interleavedInput){
        util::error_exit("<in2> is not allowed when <in1> is specified as interleaved mode by (--interleaved_in)");
    }
    // check overwrite
    if(!this->out1.empty()){
        if(this->out1 == this->out2){
            util::error_exit("read1 output (--out1) and read1 output (--out2) should be different");
        }
        if(this->donotOverwrite && util::isfile(this->out1)){
            util::error_exit(this->out1 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(!this->out2.empty()){
        if(this->donotOverwrite && util::isfile(this->out2)){
            util::error_exit(this->out2 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(this->donotOverwrite){
        if(util::isfile(this->jsonFile)){
            util::error_exit(this->jsonFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(util::isfile(this->htmlFile)){
            util::error_exit(this->htmlFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    // check compression level
    if(this->compression < 1 || this->compression > 9){
        util::error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");
    }
    // check readsToProcess 
    if(this->readsToProcess < 0){
        util::error_exit("the number of reads to process (--reads_to_process) cannot be negative");
    }
    // check thread, maximum needed is 16
    if(this->thread < 1){
        this->thread = 1;
    }else if(this->thread > 16){
        std::cerr <<  "16 threads is enough, although you specified " << thread << "\n";
        this->thread = 16;
    }
    // check trim args
    if(this->trim.front1 < 0 || this->trim.front1 > 30){
        util::error_exit("-front1 (--front1) should be 0 ~ 30, suggest 0 ~ 4");
    }
    if(this->trim.tail1 < 0 || this->trim.tail1 > 100){
        util::error_exit("-tail1 (--tail1) should be 0 ~ 100, suggest 0 ~ 4");
    }
    if(this->trim.front2 < 0 || this->trim.front2 > 30){
        util::error_exit("-front2 (--front2){ should be 0 ~ 30, suggest 0 ~ 4");
    }
    if(this->trim.tail2 < 0 || this->trim.tail2 > 100){
        util::error_exit("-tail2 (--tail2){ should be 0 ~ 100, suggest 0 ~ 4");
    }
    // to do ...//
    return true;
}

bool Options::shallDetectAdapter(bool isR2){
    if(!this->adapter.enableTriming){
        return false;
    }
    if(isR2){
        return this->isPaired() && this->adapter.enableDetectForPE && this->adapter.inputAdapterSeqR2 == "auto";
    }else{
        if(this->isPaired()){
            return this->adapter.enableDetectForPE && this->adapter.inputAdapterSeqR1 == "auto";
        }else{
            return this->adapter.inputAdapterSeqR1 == "auto";
        }
    }
}

void Options::initIndexFilter(const std::string& blacklistFile1, const std::string& blacklistFile2, int threshold){
    if(blacklistFile1.empty() && blacklistFile2.empty()){
        return;
    }
    if(!blacklistFile1.empty()){
        util::valid_file(blacklistFile1);
        this->indexFilter.blacklist1 = this->makeListFromFileByLine(blacklistFile1);
    }
    if(!blacklistFile2.empty()){
        util::valid_file(blacklistFile2);
        this->indexFilter.blacklist2 = this->makeListFromFileByLine(blacklistFile2);
    }
    if(this->indexFilter.blacklist1.empty() && this->indexFilter.blacklist2.empty()){
        return;
    }
    this->indexFilter.enabled = true;
    this->indexFilter.threshold = threshold;
}

std::vector<std::string> Options::makeListFromFileByLine(const std::string& filename){
    std::vector<std::string> ret;
    std::ifstream fr(filename);
    std::string line;
    while(std::getline(fr, line)){
        util::strip(line);
        if(line.find_first_not_of("ATCG") != std::string::npos){
            util::error_exit("processing " + filename + ", each line should be one index, which can only contain A/T/C/G");
        }
        ret.push_back(line);
    }
    return ret;
}

std::string Options::getAdapter1(){
    if(this->adapter.inputAdapterSeqR1 == "" || this->adapter.inputAdapterSeqR1 == "auto"){
        return "unspecified";
    }else{
        return this->adapter.inputAdapterSeqR1;
    }
}

std::string Options::getAdapter2(){
    if(this->adapter.inputAdapterSeqR2 == "" || this->adapter.inputAdapterSeqR2 == "auto"){
        return "unspecified";
    }else{
        return this->adapter.inputAdapterSeqR2;
    }
}

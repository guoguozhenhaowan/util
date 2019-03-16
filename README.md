utilities

|Materials|Information
|---------|-----------
|[ascii.characters.for.table](./ascii.characters.for.table)|ascii code for table drawing
|[ATCG.Binary](./ATCG.Binary)|Nucleotides to binary bits
|[plotly.js](./plotly.js)|javascripts to display motive tables
|[Doxygen](./Doxygen)|doxygen configure file for this repo

|Scripts|Functions
|-------|---------
|[getTrs.py](./getTrs.py)|get transcript fasta from genome
|[getCDSBed.py](./getCDSBed.py)|get CDS bed region in genome
|[getGenBank.py](./getGenBank.py)|get genbank file of each transcript
|[queryLib.py](./queryLib.py)|fq library cruder|  

|Sources|Contents
|----|-----------
|[tbit.cpp](./tbit.cpp)|convert any numbet to binary bits
|[extractfa.cpp](./extractfa.cpp)|extract fasta by fixed pattern in name
|[fa2bed.cpp](./fa2bed.cpp)|fasta to bed
|[fuzzy.h](./fuzzy.h)|c++ bitap search template
|[verpair.c](./verpair.c)|get different read name of two fastq file
|[versqual.c](./versqual.c)|check fq quality && sequence length
|[CLI.hpp](./CLI.hpp)|recode from https://github.com/CLIUtils/CLI11
|[dirutil.h](./dirutil.h)|recode from https://github.com/gpakosz/whereami
|[util.h](./util.h)|useful functions
|[gtutil.cpp](./gtutil.cpp)|google test for [util.h](./util.h)
|[commoh.h](./common.h)|common const definition header file
|[seq.h](./seq.h)|Seq class to represent a nucleotide sequence| 
|[gtseq.cpp](./gtseq.cpp)|google test for [seq.h](./seq.h)
|[read.h](./read.h)|Read and ReadPair class to represent fastq record
|[read.cpp](./read.cpp)|implementation of some functions in [read.h](./read.h)
|[gtread.cpp](./gtread.cpp)|google test for [read.h](./read.h)
|[fqreader.h](./fqreader.h)|FqReader and FqReaderPair class to read fastq/pair
|[fqreader.cpp](./fqreader.cpp)|implementation of some functions in [fqreader.h](./fqreader.h)
|[gtfqreader.cpp](./gtfqreader.cpp)|google test for [fqreader.h](./fqreader.h)
|[nucleotidetree.h](./nucleotidetree.h)|NucleotideNode and NucleotideTree class to represent nucleotide and nucleotide sequences
|[nucleotidetree.cpp](./nucleotidetree.cpp)|implementation of some functions in [nucleotidetree.h](./nucleotidetree.h)
|[gtnucleotidetree.cpp](./gtnucleotidetree.cpp)|google test for [nucleotidetree.h](./nucleotidetree.h)
|[knownadapters.h](./knownadapters.h)|known illumina adapters(with 3' end 1 base off)
|[evaluator.h](./evaluator.h)|Evaluator class to evaluate representation sequence, read length, read number, adapter sequence, etc|
|[evaluator.cpp](./evaluator.cpp)|implementation of some functions in [evaluator.h](./evaluator.h)
|[gtevaluator.cpp](./gtevaluator.cpp)|google test for [evaluator.h](./evaluator.h)
|[stats.h](./stats.h)|Stats class to do fastq QC, Kmer, ORSA read by read or merge Stats\* list etc
|[stats.cpp](./stats.cpp)|implementation of some functions in [stats.h](./stats.h)
|[gtstats.cpp](./gtstats.cpp)|google test for [stats.h](./stats.h)
|[overlapanalysis.h](./overlapanalysis.h)|OverlapAnalysis and OverlapResult Class to do overlap analysis and store results
|[overlapAnalysis.cpp](./overlapanalysis.cpp)|implementation of some functions in [overlapanalysis.h](./overlapanalysis.h)
|[gtoverlapanalysis.cpp](./gtoverlapanalysis.cpp)|google test for [overlapanalysis.h](./overlapanalysis.h)
|[filter.h](./filter.h)|Filter class to do various filter
|[filter.cpp](./filter.cpp)|implementation of some functions in [filter.h](./filter.h)
|[gtfilter.cpp](./gtfilter.cpp)|google test for [filter.h](./filter.h)
|[writer.h](./writer.h)|Writer class to do plain or gz format writting
|[writer.cpp](./writer.cpp)|implementation of some functions in [writer.h](./writer.h)
|[gtwriter.cpp](./gtwriter.cpp)|google test for [writer.h](./writer.h)

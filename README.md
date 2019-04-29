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
|[getGenBank.py](./getGenBank.py)|get genbank file of each transcript
|[queryLib.py](./queryLib.py)|fq library cruder|  

|Sources|Contents
|----|-----------
|[tobit.cpp](./tbit.cpp)|convert any numbet to binary bits
|[extractfa.cpp](./extractfa.cpp)|extract fasta by fixed pattern in name
|[fa2bed.cpp](./fa2bed.cpp)|fasta to bed
|[fuzzy.h](./fuzzy.h)|c++ bitap search template
|[verpair.c](./verpair.c)|get different read name of two fastq file
|[versqual.c](./versqual.c)|check fq quality && sequence length
|[CLI.hpp](./CLI.hpp)|recode from https://github.com/CLIUtils/CLI11
|[dirutil.h](./dirutil.h)|recode from https://github.com/gpakosz/whereami
|[htmlutil.h](./htmlutil.h)|html writting utilities
|[jsonutil.h](./jsonutil.h)|json writting utilities
|[util.h](./util.h)|useful functions
|[flagdel.cpp](./flagdel.cpp)|unmask some flag in a bam
|[filereader.h](./filereader.h)|general plain/gz format file reader
|[filereader.cpp](./filereader.cpp)|implementation of some functions in [filereader.h](./filereader.h)
|[filewriter.h](./filewriter.h)|general plain/gz format file writer
|[filewriter.cpp](./filewriter.cpp)|implementation of some functions in [filewriter.h](./filewriter.h)
|[getReadPairByFlag.cpp](./getReadPairByFlag.cpp)|get/filter readpair from bam alignment flag
|[cleanFaName.cpp](./cleanFaName.cpp)|clean fasta file names by remove all contens after the first blank space
|[regionStepGCDepth.gcc](./regionStepGCDepth.gcc)|step-wise gc and depth table generator for bam alignment
|[unalignedseq.h](./unalignedseq.h)|unaligned sequence class 
|[onlinebwa.h](./onlinebwa.h)|online bwa mem
|[onlinebwa.cpp](./onlinebwa.cpp)|some functions implementation of [onlinebwa.h](./onlinebwa.h)
|[extractBamByZF.cpp](./extractBamByZF.cpp)|extract bam records by ZF flag to bam files
|[splitBamByZF.cpp](./splitBamByZF.cpp)|split bam according to ZF flag to bam files
|[parseclw.cpp](./parseclw.cpp)|muscle msa result parser
|[ parsephy.cpp](./parsephy.cpp)|parse newick tree to get merged groups

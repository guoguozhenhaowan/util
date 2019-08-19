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
|[tobit.cpp](./tobit.cpp)|convert any numbet to binary bits
|[extractfa.cpp](./extractfa.cpp)|extract fasta by fixed pattern in name
|[fa2bed.cpp](./fa2bed.cpp)|fasta to bed
|[fuzzy.h](./fuzzy.h)|c++ bitap search template
|[verpair.c](./verpair.c)|get different read name of two fastq file
|[versqual.c](./versqual.c)|check fq quality && sequence length
|[CLI.hpp](./CLI.hpp)|recode from [CLI11](https://github.com/CLIUtils/CLI11)
|[json.hpp](./json.hpp)|recode from [json](https://github.com/nlohmann/json)
|[ctml.hpp](./ctml.hpp)|clone from [CTML](https://github.com/tinfoilboy/CTML)
|[dirutil.h](./dirutil.h)|recode from [whereami](https://github.com/gpakosz/whereami)
|[htmlutil.h](./htmlutil.h)|html writting utilities
|[jsonutil.h](./jsonutil.h)|json writting utilities
|[util.h](./util.h)|useful functions for common usage
|[bamutil.h](./bamutil.h)|useful functions to work with bam
|[flagdel.cpp](./flagdel.cpp)|unmask some flag in a bam
|[filereader.h](./filereader.h)|general plain/gz format file reader
|[filewriter.h](./filewriter.h)|general plain/gz format file writer
|[getReadPairByFlag.cpp](./getReadPairByFlag.cpp)|get/filter readpair from bam alignment flag
|[getAlnByRead.cpp](./getAlnByRead.cpp)|get alignment record of one read by read name
|[cleanFaName.cpp](./cleanFaName.cpp)|clean fasta file names by remove all contens after the first blank space
|[regionStepGCDepth.cpp](./regionStepGCDepth/regionStepGCDepth.cpp)|step-wise gc and depth table generator for bam alignment
|[unalignedseq.h](./unalignedseq.h)|unaligned sequence class 
|[onlinebwa.h](./onlinebwa.h)|online bwa mem
|[onlinebwa.cpp](./onlinebwa.cpp)|some functions implementation of [onlinebwa.h](./onlinebwa.h)
|[extractBamByZF.cpp](./extractBamByZF.cpp)|extract bam records by ZF flag to bam files
|[splitBamByZF.cpp](./splitBamByZF.cpp)|split bam according to ZF flag to bam files
|[parseclw.cpp](./parseclw.cpp)|muscle msa result parser
|[parsephy.cpp](./parsephy.cpp)|parse newick tree to get merged groups
|[kmerStat.cpp](./kmerStat.cpp)|count kmers in a fasta/fastq file
|[kseq.h](./kseq.h)|customized version kseq.h which seperate read and store
|[compareBamRead.cpp](./compareBamRead.cpp)|compare reads included in two bams
|[getContig.cpp](./getContig.cpp)|get region of a contig in reference
|[statutil.h](./statutil.h)|utilities to calculate some statistical items
|[intervalTree.h](./intervalTree.h)|a simple interval tree support construct once and query later
|[samheader.h](./samheader.h)|utilities to operate on sam/bam header, revised from[samtools](https://github.com/samtools/samtools)
|[samheader.cpp](./samheader.cpp)|some functions inmplmentation of [samheader.h](./samheader.h), revised from[samtools](https://github.com/samtools/samtools)
|[simulateSV.cpp](./simulateSV.cpp)|simulate SV from reference
|[fq2fa.c](./fq2fa.c)|a simple tool to convert fastq to fasta
|[modifyGTF.cpp](./modifyGTF.cpp)|modify UCSC gtf file to add gene name in group column
|[getFeatureTsv.cpp](./getFeatureTsv.cpp)|extract feature tsv file from UCSC ref gene tsv file with accession number
|[getSUR.c](./getSUR.c)|extract single unmapped bam records
|[getAlnByZF.cpp](./getAlnByZF.cpp)|extract bam record by ZF tag
|[matchFqSeq.c](./matchFqSeq.c)|extract read sequence exactly containning provided seq

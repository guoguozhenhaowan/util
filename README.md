utilities

|Materials|Information
|---------|-----------
|[ascii.characters.for.table](./ascii.characters.for.table)|ascii code for table drawing
|[ATCG.Binary](./ATCG.Binary)|Nucleotides to binary bits
|[getTrs.py](./getTrs.py)|get transcript fasta from genome
|[queryLib.py](./queryLib.py)|fq library cruder

|Scripts|Functions
|-------|---------
|[getCDSBed.py](./getCDSBed.py)|get CDS bed region in genome
|[getGenBank.py](./getGenBank.py)|get genbank file of each transcript
|[queryLib.py](./queryLib.py)|fq library cruder|  

|Sources|Contents
|----|-----------
|[extractfa.cpp](./extractfa.cpp)|extract fasta by fixed pattern in name
|[fa2bed.cpp](./fa2bed.cpp)|fasta to bed
|[fuzzy.h](./fuzzy.h)|c++ bitap search template
|[verpair.c](./verpair.c)|get different read name of two fastq file
|[versqual.c](./versqual.c)|check fq quality && sequence length
|[CLI.hpp](./CLI.hpp)|recode from https://github.com/CLIUtils/CLI11
|[dirutil.h](./dirutil.h)|recode from https://github.com/gpakosz/whereami
|[util.h](./util.h)|useful functions
|[gtutil.cpp](./gtutil.cpp)|google test for [util.h](./util.h)
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


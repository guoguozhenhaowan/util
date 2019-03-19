# test util.h 
g++ -std=c++11 -lgtest gtutil.cpp
# test seq.h
g++ -std=c++11 -lgtest gtseq.cpp
# test read.h
g++ -std=c++11 -lgtest gtread.cpp read.cpp
# test fqreader.h
g++ -std=c++11 -lgtest -lz gtfqreader.cpp fqreader.cpp
# test polyx.h
g++ -std=c++11 -lgtest gtpolyx.cpp polyx.cpp
# test writer.h
g++ -std=c++11 -lz -lgtest gtwriter.cpp writer.cpp fqreader.cpp
# test nucleotidetree.h
g++ -std=c++11 -lgtest gtnucleotidetree.cpp nucleotidetree.cpp
# test overlapanalysis.h
g++ -std=c++11 -lgtest gtoverlapanalysis.cpp overlapanalysis.cpp
# test filter.h
g++ -std=c++11 -lgtest gtfilter.cpp filter.cpp read.cpp options.cpp
# test duplicate.h
g++ -std=c++11 -lgtest -lz gtduplicate.cpp duplicate.cpp fqreader.cpp options.cpp
# test basecorrector.h
g++ -std=c++11 -lgtest gtbasecorrector.cpp basecorrector.cpp filterresult.cpp overlapanalysis.cpp
# test adaptertrimmer.h
g++ -std=c++11 -lgtest gtadaptertrimmer.cpp adaptertrimmer.cpp filterresult.cpp overlapanalysis.cpp
# test evaluator.h
g++ -std=c++11 -lgtest -lz gtevaluator.cpp evaluator.cpp read.cpp fqreader.cpp nucleotidetree.cpp options.cpp
# test stats.h
g++ -std=c++11 -lgtest -lz gtstats.cpp read.cpp stats.cpp options.cpp evaluator.cpp nucleotidetree.cpp fqreader.cpp
# test filterresult.h
g++ -std=c++11 -lgtest -lz gtfilterresult.cpp filterresult.cpp stats.cpp options.cpp evaluator.cpp nucleotidetree.cpp fqreader.cpp

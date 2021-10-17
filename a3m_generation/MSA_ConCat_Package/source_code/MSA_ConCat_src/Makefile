

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all:
	$(CXX) A3M_NoGap.cpp $(FLAG) -o A3M_NoGap
	$(CXX) A3M_SpecBloc.cpp $(FLAG) -o A3M_SpecBloc
	$(CXX) MSA_ConCat.cpp $(FLAG) -o MSA_ConCat
	cp A3M_NoGap A3M_SpecBloc ../../bin/
	cp MSA_ConCat ../../
	mv A3M_NoGap A3M_SpecBloc MSA_ConCat ../
# clean all the temp files
clean:


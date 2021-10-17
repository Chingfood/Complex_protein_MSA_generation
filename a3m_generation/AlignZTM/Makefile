SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all: align calcMeff

align: align.cpp
	$(CXX) align.cpp $(FLAG) -o align

calcMeff: calcMeff.cpp
	$(CXX) calcMeff.cpp $(FLAG) -o calcMeff

clean:
	rm align calcMeff
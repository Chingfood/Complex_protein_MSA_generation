

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all: tgt cat util

tgt:
	# --- A3M_To_TGT ---- #
	cd MSA_To_PSSM_cdhit_src&&make
	mv MSA_To_PSSM_cdhit_src/MSA_To_PSSM ../util
	$(CXX) A3M_To_TGT.cpp $(FLAG) -o A3M_To_TGT
	mv A3M_To_TGT ../

cat:
	# --- MSA_ConCat ---- #
	cd MSA_ConCat_src&&make
	mv A3M_NoGap A3M_SpecBloc ../bin
	mv MSA_ConCat ../

util:
	# ---- util ---- #
	cd util_src&&make
	#-> for A3M_To_TGT
	mv Verify_FASTA A3M_To_A2M_seqid A3M_To_PSI SS8_To_SS3 ACC_ReFormat ../util

# clean all the temp files
clean:


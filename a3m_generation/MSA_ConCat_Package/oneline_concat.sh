#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: ./oneline_concat.sh <input_fasta_1> <input_fasta_2> [CPU_NUM] "
	exit
fi
CurRoot="$(pwd)"

# ------- CPU number ------ #
BindX_CPU=1
if [ $# -gt 2 ]
then
	BindX_CPU=$3
fi

# ------ part 0 ------ # related path
infile1=$1
infile2=$2
fulnam1=`basename $1`
relnam1=${fulnam1%.*}
fulnam2=`basename $2`
relnam2=${fulnam2%.*}

#-> data path
species_data=data/names.dmp_species_taxid_better_fixed3_TaxTree

# ------ part 1 ------ # generate A3M
tmp=tmp/$relnam1"_"$relnam2"_MSA_ConCat"
mkdir -p $tmp
echo "part 1.1 : generate A3M file for the 1st input sequence $relnam1"
./Fast_TGT.sh -i $infile1 -c $BindX_CPU -o $tmp 1> ws1 2> ws2
rm -f ws1 ws2
echo "part 1.2 : generate A3M file for the 2nd input sequence $relnam2"
./Fast_TGT.sh -i $infile2 -c $BindX_CPU -o $tmp 1> ws1 2> ws2
rm -f ws1 ws2

# ----- part 2 ------- # process A3M
echo "part 2.1 : process the 1st A3M file $relnam1.a3m"
bin/A3M_NoGap $tmp/$relnam1.a3m $tmp/$relnam1.a3m_nogap 1>ws1 2>ws2
bin/A3M_SpecBloc $tmp/$relnam1.a3m_nogap $species_data $tmp/$relnam1.a3m_specbloc 1> $tmp/$relnam1.block 2>ws2
echo "part 2.2 : process the 2nd A3M file $relnam2.a3m"
bin/A3M_NoGap $tmp/$relnam2.a3m $tmp/$relnam2.a3m_nogap 1>ws1 2>ws2
bin/A3M_SpecBloc $tmp/$relnam2.a3m_nogap $species_data $tmp/$relnam2.a3m_specbloc 1> $tmp/$relnam2.block 2>ws2
rm -f ws1 ws2

# ----- part 3 ------- # concatenate two MSAs
echo "part 3 : concatenate two MSAs"
./MSA_ConCat $tmp/$relnam1.a3m_specbloc $tmp/$relnam2.a3m_specbloc $tmp/${relnam1}_${relnam2}.a3m_concat 1>ws1 2>ws2
rm -f ws1 ws2




#!/bin/bash

if [ $# -ne 5 ]
then
	echo "Usage: ./genMSA.sh <infile> <outfolder> <coverage> <arg4hh> <Meff>"
	exit
fi
path=`pwd`
infile=`readlink -f $1`
outfolder=`readlink -f $2`
coverage=$3
arg4hh=$4
Meff=$5

TMPDIR=$outfolder

mkdir -p $outfolder
if [ ! -f $infile ]; then echo "infile $infile not found !!" >&2; exit 1; fi
if [ ! -d $outfolder ]; then echo "outfolder $outfolder not folder" >&2; exit 1; fi


fulnam=`basename $infile`
wsrelnam=${fulnam%.*}


relnam=${infile%.*}
relnam=${relnam##*/}
orelnam=$outfolder/$relnam
relnam=$TMPDIR/$relnam

fasta_raw_file=$relnam.fasta_raw
seq_file=$relnam.seq
a3m_file=$relnam.a3m
a3m_nogap_file=$relnam.a3m_nogap
a3m_specbloc_file=$relnam.a3m_specbloc

ofasta_raw_file_compress=$orelnam.fasta_raw.tar.gz
oseq_file_compress=$orelnam.seq.tar.gz
oa3m_file_compress=$orelnam.a3m.tar.gz
oa3m_specbloc_file_compress=$orelnam.a3m_specbloc.tar.gz


cp $infile $fasta_raw_file
cd $TMPDIR

function o2t { if [ ! -f $2.tar.gz ]; then cp $1 $2.tar.gz; fi; tar -xzf $2.tar.gz; }
function t2o { tar -czf $1.tar.gz ${1##*/}; if [ ! -f $2 ]; then cp $1.tar.gz $2; fi }
function cpback { if [[ -f $1 && ! -f $2 ]]; then t2o $1 $2; fi }

if [ ! $MSA_ConCat_Package ]; then echo "MSA_ConCat_Package not specified"; exit; fi
$MSA_ConCat_Package/util/Verify_FASTA $fasta_raw_file $seq_file

# ----- determine coverage ---- #
if [ $coverage -eq -1 ]; then a=60; b=`tail -n1 $seq_file | wc | awk '{print int(7000/($3-1))}'`; if [ $a -gt $b ]; then coverage=$b; else coverage=$a; fi; fi

# ---- generate A3M file -------- #
if [ ! -f $oa3m_file_compress ]; then
	if [ ! $HHSUITE ]; then echo "HHsuite not specified"; exit; fi
	echo "hhblits start with parameters -cov $coverage $arg4hh"
#	export HHLIB=$HHSUITE/lib/hh
#	export LD_LIBRARY_PATH=/opt/hhsuite/lib:$LD_LIBRARY_PATH
	$HHSUITE/bin/hhblits -i $seq_file -o $relnam.hhr -oa3m $a3m_file -cov $coverage $arg4hh  1> $relnam.hho 2> $relnam.hhe
#		-i $seq_file -o $relnam.hhr -oa3m $relnam.a3m -d /mnt/data/conmod_databases/uniprot20_series/uniprot20_2016_02/$uniprot20	\
#		-e $eval -cpu $cpu_num -n $iteration -maxfilt 500000 -diff inf -id 99 -cov $coverage -nodiff
#	rm $relnam.hhr
	echo "hhblits done"
else
	o2t $oa3m_file_compress $a3m_file
fi
cpback $a3m_file $oa3m_file_compress a3m

# delete these lines
#if [[ $NSLOTS != 4 ]]; then t2o $a3m_file $oa3m_file_compress; exit; fi


if [ ! -f $oa3m_specbloc_file_compress ]; then
	echo "process the A3M file $relnam.a3m"
	$MSA_ConCat_Package/bin/A3M_NoGap $a3m_file $a3m_nogap_file 1>ws1 2>ws2
	rm -f ws1 ws2
#	cp $MSA_ConCat_Package/data/names.dmp_species_taxid_better_fixed3_TaxTree.tar.gz $TMPDIR/
#	tar -xzf $TMPDIR/names.dmp_species_taxid_better_fixed3_TaxTree.tar.gz
#	rm -f $TMPDIR/names.dmp_species_taxid_better_fixed3_TaxTree.tar.gz
	species_data=$MSA_ConCat_Package/data/names.dmp_species_taxid_better_fixed3_TaxTree
	$MSA_ConCat_Package/bin/A3M_SpecBloc $a3m_nogap_file $species_data $a3m_specbloc_file 1> $relnam.block 2>ws2
	rm -f ws2
fi




cpback $seq_file $oseq_file_compress seq
cpback $a3m_file $oa3m_file_compress a3m
cpback $a3m_specbloc_file $oa3m_specbloc_file_compress a3m_specbloc

#!/bin/bash

function usage() { echo './align.sh <-a A3M a> <-b A3M b> <-o outrelnam> [-c coverage] [-i identity] [-m min_delta_gene] [-M max_delta_gene] [-j joint method] [-d GPUdevice] [-T TMPDIR] [-e Meff] [-t #threads]'; exit 1; }

infile1=''
infile2=''
relnam=''
min_dgene=1
max_dgene=20
GPUdevice=-2
coverage=-1
identity=90
Meff=.7
nCPU=1

while getopts ":a:b:o:c:i:m:M:j:d:T:e:t:" opt; do case $opt in
	a) infile1=`readlink -f $OPTARG`;;
	b) infile2=`readlink -f $OPTARG`;;
	o) relnam=`readlink -f $OPTARG`;;
	c) coverage=$OPTARG;;
	i) identity=$OPTARG;;
	m) min_dgene=$OPTARG;;
	M) max_dgene=$OPTARG;;
	j) joint_method=$OPTARG;;
	d) GPUdevice=$OPTARG;;
	T) TMPDIR=`readlink -f $OPTARG`;;
	e) Meff=$OPTARG;;
	t) nCPU=$OPTARG;;

	\?) echo "Invalid option: -$OPTARG" >&2; usage; exit 1;;
	:) echo "Option -$OPTARG requires an argument." >&2; usage; exit 1;;
esac; done
#echo "GPUdevice $GPUdevice"
#echo "nCPU $nCPU"



if [[ ! -f $infile1.tar.gz ]]; then echo "infile1 $infile1 not found"; exit 1; fi
if [[ ! -f $infile2.tar.gz ]]; then echo "infile2 $infile2 not found"; exit 1; fi
outfolder=${relnam%/*}
if [[ ! $TMPDIR ]]; then TMPDIR=$outfolder; fi
mkdir -p $outfolder
if [[ ! -d $outfolder ]]; then echo "outfolder $outfolder not folder"; exit 1; fi
cp $infile1.tar.gz $infile2.tar.gz $TMPDIR/


fulnam1=`basename $infile1`
wsrelnam1=${fulnam1%.*}
fulnam2=`basename $infile2`
wsrelnam2=${fulnam2%.*}
wsrelnam="${wsrelnam1}_${wsrelnam2}"

orelnam=$relnam
relnam=$TMPDIR/${relnam##*/}

cd $TMPDIR

infile1=${infile1##*/}
infile2=${infile2##*/}
tar -xzf $infile1.tar.gz
tar -xzf $infile2.tar.gz
rm -f $infile1.tar.gz $infile2.tar.gz

a3m_file_raw=$relnam.a3m_raw
a3m_file=$relnam.a3m
seq_file=$relnam.seq

oa3m_file_raw_compress=$orelnam.a3m_raw.tar.gz
oa3m_file_compress=$orelnam.a3m.tar.gz
oseq_file_compress=$orelnam.seq.tar.gz


function o2t { if [ ! -f $2.tar.gz ]; then cp $1 $2.tar.gz; fi; tar -xzf $2.tar.gz; }
function t2o { tar -czf $1.tar.gz ${1##*/}; if [ ! -f $2 ]; then cp $1.tar.gz $2; fi }

if [ ! -f $oa3m_file_raw_compress ]; then database=$AlignZTM/id2loc; $AlignZTM/align $infile1 $infile2 $min_dgene $max_dgene $database $joint_method $a3m_file_raw;
else o2t $oa3m_file_raw_compress $a3m_file_raw; fi

if [[ -f $a3m_file_raw && ! -f $oa3m_file_compress ]]; then 
	if [[ $coverage -eq -1 ]]; then a=60; b=`head -n2 $a3m_file_raw | tail -n1 | wc | awk '{print int(7000/($3-1))}'`; if [ $a -gt $b ]; then coverage=$b; else coverage=$a; fi; fi
	$HHSUITE/bin/hhfilter -i $a3m_file_raw -cov $coverage -id $identity -o $a3m_file; 
elif [[ -f $oa3m_file_compress ]]; then o2t $oa3m_file_compress $a3m_file; fi


if [[ -f $a3m_file_raw && ! -f $oa3m_file_raw_compress ]]; then t2o $a3m_file_raw $oa3m_file_raw_compress; fi
if [[ -f $a3m_file && ! -f $oa3m_file_compress ]]; then t2o $a3m_file $oa3m_file_compress; fi
if [[ -f $seq_file && ! -f $oseq_file_compress ]]; then t2o $seq_file $oseq_file_compress; fi


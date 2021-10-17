#!/bin/bash

## This script calls UniprotID2EmblID_new.py
## Multi process for using Uniprot.org cross reference service to convert UniprotID to EMBL ID. 


ID_list=$1
source=$2
contig_list_folder=$3
mkdir -p $contig_list_folder

MY_PATH="`dirname \"$0\"`"

from=$4
to=$5
num_thread=$6
n_per_thread=$(($to-$from+$num_thread-1))
n_per_thread=$(($n_per_thread/$num_thread))

for (( i=0; i<$num_thread; i++ ))
do
	l=$((i*n_per_thread+$from))
	if [[ $l -ge $to ]]; then break; fi
	r=$((i*n_per_thread+$from+n_per_thread))
	outfile=$contig_list_folder/${l}_${r}
	if [[ $source == UniProtKB/TrEMBL ]]; then outfile=${outfile}_UniProtKB_TrEMBL;
	elif [[ $source == UniProtKB/Swiss-Prot ]]; then outfile=${outfile}_UniProtKB_Swiss-Prot;
	else outfile=${outfile}_${source};
	fi
	logfile=$outfile.log
	errfile=$outfile.err
	echo python3 download_xref.py $ID_list $source $outfile $l $r
	nohup python3 $MY_PATH/util/UniprotID2EmblID_new.py $ID_list $source $outfile $l $r > $logfile 2>$errfile  &
	sleep 3
done
wait

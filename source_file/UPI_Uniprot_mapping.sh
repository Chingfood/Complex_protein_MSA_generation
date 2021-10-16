#!/bin/bash
## This bash script call UniRefID2UniProtID_parallel.py
## Converts UnirefID to UniprotID using multiple processes

ID_list=$1
uniprot_list_folder=$2
mkdir -p $uniprot_list_folder

MY_PATH="`dirname \"$0\"`"

from=$3
to=$4
num_thread=$5
n_per_thread=$(($to-$from+$num_thread-1))
n_per_thread=$(($n_per_thread/$num_thread))

for (( i=0; i<$num_thread; i++ ))
do
	l=$((i*n_per_thread+$from))
	if [[ $l -ge $to ]]; then break; fi
	r=$((i*n_per_thread+$from+n_per_thread))
	outfile=$uniprot_list_folder/${l}_${r}.out
	logfile=$uniprot_list_folder/${l}_${r}.log
	errfile=$uniprot_list_folder/${l}_${r}.err
	echo python3 UniRefID_UniProtID_mapping_parallel.py $ID_list $outfile $l $r
	nohup python3 $MY_PATH/util/UniRefID_UniProtID_mapping_parallel.py $ID_list $outfile $l $r > $logfile 2>$errfile  &
	sleep 3
done
wait

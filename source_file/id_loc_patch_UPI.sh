#!/bin/bash
## This bash script call UniRefID2UniProtID_parallel.py
## Converts UnirefID to UniprotID using multiple processes

id_loc=$1
UPI_UNI_mapping=$2
id_loc_patch_folder=$3
mkdir -p $id_loc_patch_folder

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
	outfile=$id_loc_patch_folder/${l}_${r}.out
	logfile=$id_loc_patch_folder/${l}_${r}.log
	errfile=$id_loc_patch_folder/${l}_${r}.err
	echo python3 id_loc_patch_UPI.py $id_loc $UPI_UNI_mapping $outfile $l $r
	nohup python3 $MY_PATH/util/id_loc_patch_UPI_parallel.py $id_loc $UPI_UNI_mapping $outfile $l $r > $logfile 2>$errfile  &
	sleep 3
done
wait

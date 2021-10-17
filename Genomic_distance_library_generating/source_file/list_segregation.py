#!/bin/bash
## This bash script call UniRefID2UniProtID_parallel.py
## Converts UnirefID to UniprotID using multiple processes

ID_list=$1

MY_PATH="`dirname \"$0\"`"

from=$2
to=$3
num_thread=$4
n_per_thread=$(($to-$from+$num_thread-1))
n_per_thread=$(($n_per_thread/$num_thread))

for (( i=0; i<$num_thread; i++ ))
do
	l=$((i*n_per_thread+$from))
	if [[ $l -ge $to ]]; then break; fi
	r=$((i*n_per_thread+$from+n_per_thread))
	echo "l is $l"
	echo "r is $r"
	sleep 3
done

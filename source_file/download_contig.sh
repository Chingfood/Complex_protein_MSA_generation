#!/bin/bash

ID_list=$1
outfolder=$2

from=$3
to=$4
num_thread=$5
n_per_thread=$(($to-$from+$num_thread-1))
n_per_thread=$(($n_per_thread/$num_thread))

MY_PATH="`dirname \"$0\"`"

for (( i=0; i<$num_thread; i++ ))
do
	l=$((i*n_per_thread+$from))
	if [[ $l -ge $to ]]; then break; fi
	r=$((i*n_per_thread+$from+n_per_thread))
	outfile=$outfolder${l}_${r}
	logfile=$outfile.log
	echo python3 download_contig.py $ID_list $outfile $l $r
#	exit
	nohup python3 $MY_PATH/util/download_contig.py $ID_list $outfile $l $r > $logfile &
#	exit
	sleep 3
done
wait

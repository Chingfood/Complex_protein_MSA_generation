#!/bin/bash

# ./run.sh uniclust30_2017_10_a3m.ffdata 2017_10
#/mnt/data/conmod_databases/uniclust30_series/uniclust30_2018_08/uniclust30_2018_08_a3m.ffdata

#use python 3.8
#The following code activate python 3.8. If you do not require this, comment out the following two lines.

source /home/qingyliu/miniconda3/etc/profile.d/conda.sh
conda activate python3

a3m_ffindex=$1
version=$2

MY_PATH="`dirname \"$0\"`"

mkdir -p $MY_PATH/temp/data/$version

#a3m_ffindex=uniclust30_${version}_a3m.ffdata
uniref_ID_list=$MY_PATH/temp/data/$version/uniref_ID_list_${version}
uniprot_ID_list=$MY_PATH/temp/data/$version/uniprot_ID_list_${version}
xref_list_folder=$MY_PATH/temp/data/$version/xref_list_${version}
xref_batch_folder=$MY_PATH/temp/data/$version/xref_batch_${version}
xref_folder=$MY_PATH/temp/data/$version/xref_${version}
merged_xref_list=$MY_PATH/temp/data/$version/merged_xref_list_${version}
contig_list=$MY_PATH/temp/data/$version/contig_list_${version}
wgs_list=$MY_PATH/temp/data/$version/wgs_list_${version}
CDS_folder=$MY_PATH/temp/data/$version/CDS_folder_${version}
CDS_wgs_folder=$MY_PATH/temp/data/$version/CDS_wgs_folder_${version}
loc_folder=$MY_PATH/temp/data/$version/loc_folder_${version}
id_loc=$MY_PATH/temp/data/$version/id_loc_${version}
mkdir -p $xref_list_folder $xref_batch_folder $xref_folder $CDS_folder $CDS_wgs_folder $loc_folder

# extract IDs in our database
                                   
#echo "Collecting UniRef IDs"
#python3 $MY_PATH/source_file/util/extract_ID_UniRef.py $a3m_ffindex $uniref_ID_list

#python3 $MY_PATH/source_file/util/UniRefID_UniProtID_seperate.py $uniref_ID_list 

##May write it to have multiple process to speed up this step.
#num_query=`wc -l $uniref_ID_list | cut -d' ' -f1`
#python3 $MY_PATH/source_file/util/UniRefID2UniProtID.py $uniref_ID_list 

#cat $MY_PATH/temp/data/$version/uniref_id_UNI $MY_PATH/temp/data/$version/uniref_id_UPI_converted > $MY_PATH/temp/data/$version/uniprot_ID_list_redundant
#python3 $MY_PATH/source_file/util/file_line_redundancy_remove.py $MY_PATH/temp/data/$version/uniprot_ID_list_redundant $uniprot_ID_list

#rm $MY_PATH/temp/data/$version/uniref_id_UNI 
#rm $MY_PATH/temp/data/$version/uniref_id_UPI
#rm $MY_PATH/temp/data/$version/uniref_id_UPI_converted
#rm $MY_PATH/temp/data/$version/uniprot_ID_list_redundant


#num_query=`wc -l $uniprot_ID_list | cut -d' ' -f1`
#echo "We have $num_query queries in total"
## 97 884 147 in total
#source=UniProtKB/TrEMBL
## IF NOT IN TrEMBL, swissprot is checked as well in download_xref.py
#echo "Downloading $source entries"
#from=0; to=$num_query;  num_thread=30   # left-close, right-open
##from=0;        to=10000;       num_thread=20   # left-close, right-open
##i=0;   k=7000000;      from=$(($i*$k)); to=$(($(($i+1))*$k))
#./$MY_PATH/source_file/download_xref.sh $uniprot_ID_list $source $xref_list_folder $from $to $num_thread
#
#mkdir -p $xref_list_folder/temp
#mv $xref_list_folder/*.log $xref_list_folder/temp 2>/dev/null
#mv $xref_list_folder/*.err $xref_list_folder/tem1p 2>/dev/null
#cat $xref_list_folder/*EMBL >  $merged_xref_list
#

## To verify data class
#echo "verifying data"
#cut -f 5 $merged_xref_list > $MY_PATH/temp/data/$version/tmp_contig
#cut -f 4 $merged_xref_list > $MY_PATH/temp/data/$version/tmp_coding
#wc -l $MY_PATH/temp/data/$version/tmp_contig $MY_PATH/temp/data/$version/tmp_coding

########redundant check, no need####################
######python3 -c "from enaBrowserTools.python3.utils import is_wgs_sequence, is_sequence
######with open('tmp_contig', 'r') as f: assert (is_sequence(_) or is_wgs_sequence(_) for _ in [l[:-1] for l in set(f)])" &
######
######python3 -c "from enaBrowserTools.python3.utils import is_coding
######with open('tmp_coding', 'r') as f: assert all(is_coding(_) for _ in f.read().strip().split())" &
######wait
######echo "verified"
#############end####################

# To separate contig and wgs
#echo "separating contig and wgs"
#python3 $MY_PATH/source_file/util/extract_contig_wgs_ID.py $MY_PATH/temp/data/$version/tmp_contig $contig_list $wgs_list
#wc -l $contig_list $wgs_list


#rm -f $MY_PATH/temp/data/$version/tmp_contig $MY_PATH/temp/data/$version/tmp_coding
#ID_list=$contig_list;   outfolder=$CDS_folder/contig_
#num_query=`wc -l $ID_list | cut -d' ' -f1`
#echo "We have $num_query queries in total"
## 9 105 361
#from=0; to=$num_query;  num_thread=100   # left-close, right-open
#./$MY_PATH/source_file/download_contig.sh $ID_list $outfolder $from $to $num_thread
#wait
#
ID_list=$wgs_list;              outfolder=$CDS_wgs_folder/wgs_
num_query=`wc -l $ID_list | cut -d' ' -f1`
echo "We have $num_query queries in total"
# 16 704
from=0; to=$num_query;  num_thread=100   # left-close, right-open
./$MY_PATH/source_file/download_contig.sh $ID_list $outfolder $from $to $num_thread
wait
#
#rm -r $CDS_folder/temp
#rm -r $CDS_wgs_folder/temp
#mv $CDS_wgs_folder/* $CDS_folder
#
#echo "extracting CDS from raw data"
#cut -f 2,4 $merged_xref_list > tmp_merged_xref_list
#python3 $MY_PATH/source_file/util/extract_CDS_parallel.py $MY_PATH/temp/data/$version/tmp_merged_xref_list $CDS_folder $loc_folder 2> $MY_PATH/temp/data/$version/extract_CDS_$version.log
#rm tmp_merged_xref_list
#
#echo "generating id_loc file"
##python2 genFile.py.bak $loc_folder $id_loc
#python3 $MY_PATH/source_file/util/genFile.py $loc_folder $id_loc
#
#

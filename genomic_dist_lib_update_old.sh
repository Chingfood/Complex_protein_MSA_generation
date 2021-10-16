#!/bin/bash

# ./run.sh -t temp UniRef30_2020_06_a3m.ffdata 2020_06

#use python 3.8 or later and use numpy
#use enaBrowserTools. Please have the update the enaBrowsertools to the most up-to-date version.

############################# Step 1: Initiation ####################################

# the directory of the script
MY_PATH="`dirname \"$0\"`"

# the directory where it is executed
DIR=`pwd`
# the temp directory used, within $DIR
# omit the -p parameter to create a temporal directory in the default location           
# remove it first in case user provided a temporary directory
TEMP_DIR=`mktemp -d -p "$DIR"`
rm -r $TEMP_DIR


### May put enaBrowsertools path here in getopt as well
function Usage 
{
	echo $0 "[ -t TempDir ] Uniref.ffdata Uniref_time_version "
	echo "	This script updates the library for calculating genomic distance for MSA."
	echo "	Uniref.ffdata: The sequence database file. Uniref database is used for HHblits. All the protein ID from the database will be extracted. Normally it is under ffdata extension. An example be: UniRef30_2020_06_a3m.ffdata." 
	echo "	Uniref_time_version: The version of the Uniref library. If the library is updated at June 2020, then the Uniref_time_version would be 2020_06."
	echo "	TempDir: A temporary directory that store all the intermediate output and final output. If not specified, a random temporary directory will be generated automatically. "

}

while getopts ":t:" opt; do
        case ${opt} in
                t )
                  TEMP_DIR=$OPTARG
                  ;;
                \? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  Usage
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  Usage
                  exit 1
                  ;;
        esac
done

shift $((OPTIND -1))

if [ $# -ne 2 ]; then
	echo "ERROR: Not all input specified"
        Usage
	exit 1
fi

a3m_ffindex=$1
if [ ! -f $a3m_ffindex ]; then
	echo "ERROR: invalid protein library file $a3m_ffindex"
        Usage
	exit 1
fi

version=$2

mkdir -p $TEMP_DIR
echo "Temprory Directory: $TEMP_DIR"
WORK_DIR=$TEMP_DIR/$version

# a3m_ffindex=uniclust30_${version}_a3m.ffdata
uniref_ID_list=$WORK_DIR/uniref_ID_list_${version}
uniprot_list_folder=$WORK_DIR/uniprot_list_${version}
uniprot_ID_list=$WORK_DIR/uniprot_ID_list_${version}
xref_list_folder=$WORK_DIR/xref_list_${version}
merged_xref_list=$WORK_DIR/merged_xref_list_${version}
contig_list=$WORK_DIR/contig_list_${version}
wgs_list=$WORK_DIR/wgs_list_${version}
CDS_folder=$WORK_DIR/CDS_folder_${version}
CDS_wgs_folder=$WORK_DIR/CDS_wgs_folder_${version}
loc_folder=$WORK_DIR/loc_folder_${version}
id_loc=$WORK_DIR/id_loc_${version}
mkdir -p $uniprot_list_folder $xref_list_folder $CDS_folder $CDS_wgs_folder $loc_folder

################# Step 2: extract protein ID from Uniref library ###########

## extract IDs in our database
#                                   
#echo "Collecting UniRef IDs"
#python3 $MY_PATH/source_file/util/extract_ID_UniRef.py $a3m_ffindex $uniref_ID_list
#
## seperate out unirefID with regular uniprotID and unirefID that starts with UPI.
## for the ID format of the UnirefID, please visit: https://www.uniprot.org/help/uniref
#python3 $MY_PATH/source_file/util/UniRefID_UniProtID_seperate.py $uniref_ID_list $WORK_DIR/uniref_id_UPI $WORK_DIR/uniref_id_UNI  
awk '{printf ("%s\t%s\n", $1, $1) }' $WORK_DIR/uniref_id_UNI > $WORK_DIR/uniref_id_UNI_UNI


# Multi_process_version to convert UnirefID that starts with UPI to UniprotID
# To check on the process finishing situation, helper script: script/last_line_log.py may help.
#num_query=`wc -l $WORK_DIR/uniref_id_UPI | cut -d' ' -f1`
#echo $num_query
#from=0; to=$num_query;  num_thread=20   # left-close, right-open
#bash $MY_PATH/source_file/UPI_Uniprot_mapping.sh $WORK_DIR/uniref_id_UPI $uniprot_list_folder $from $to $num_thread
#wait

mkdir -p $uniprot_list_folder/temp
mv $uniprot_list_folder/*.log $uniprot_list_folder/temp 
mv $uniprot_list_folder/*.err $uniprot_list_folder/temp 
cat $uniprot_list_folder/*.out > $WORK_DIR/uniref_id_UPI_UNI_multiprocess


#cat $WORK_DIR/uniref_id_UNI_UNI $WORK_DIR/uniref_id_UPI_UNI_multiprocess > $WORK_DIR/uniref_id_UPI_UNI_full
#awk '{print $2}' $WORK_DIR/uniref_id_UPI_UNI_full > $WORK_DIR/uniprot_ID_list_redundant
#python3 $MY_PATH/source_file/util/file_line_redundancy_remove.py $WORK_DIR/uniprot_ID_list_redundant $uniprot_ID_list
##
#
#
#
################ Step 3: Cross reference to find EMBL ID from UniprotID ############
#
#
## Use the cross reference service to find the EMBL_ID of each UniprotID
## There are two methods available to do this task. 
## First method use the cross reference service maintained by ENA to find EMBL ID
## Second method use the cross reference service maintained by Uniprot.org
## Both methods are implemented here in case one not working.
## If the downloading was interupted and did not finish,  script/xref_restart.py may help with restart the downloads from where it ended. It is still recommended to re-start the download from the scratch.
## To check on the process finishing situation, helper script: script/last_line_log.py may help see if the final index is the range. This does not apply to the last chunck..
#
#num_query=`wc -l $uniprot_ID_list | cut -d' ' -f1`
#echo "We have $num_query queries in total"
## 97 884 147 in total
#source=UniProtKB/TrEMBL
## IF NOT IN TrEMBL, swissprot is checked as well in download_xref.py
#echo "Downloading $source entries"
#from=0; to=$num_query;  num_thread=50   # left-close, right-open
##from=0;        to=10000;       num_thread=20   # left-close, right-open
#bash $MY_PATH/source_file/download_xref_method1.sh $uniprot_ID_list $source $xref_list_folder $from $to $num_thread
#wait
#
#### thie methods produce a different format with only two columns, uniprot id and embl id
####bash $MY_PATH/source_file/download_xref_method2.sh $uniprot_ID_list $source $xref_list_folder $from $to $num_thread
#wait
#
#mkdir -p $xref_list_folder/temp
#mv $xref_list_folder/*.log $xref_list_folder/temp 2>/dev/null
#mv $xref_list_folder/*.err $xref_list_folder/temp 2>/dev/null
#cat $xref_list_folder/*EMBL >  $merged_xref_list
#
#
## To verify data class
#echo "verifying data"
#cut -f 5 $merged_xref_list > $WORK_DIR/tmp_contig
#cut -f 4 $merged_xref_list > $WORK_DIR/tmp_coding
#wc -l $WORK_DIR/tmp_contig $WORK_DIR/tmp_coding
#
#
################ Step 4: Downloading contig file #########################
#
#
##
## The following steps and the step after require a package called enaBrowserTools
## This package can be downloaded from: https://github.com/enasequence/enaBrowserTools or https://github.com/enasequence/enaBrowserTools.git
## Make sure: 1. This tool is up-to-date. 2. If you planning to run this step at server. Have this tool unzipped at the server. 3. Change its uncompressed folder's name to enaBrowserTools. (The default is enaBrowserTools-master) 4. copy this folder to: source_file/util/
## If the process got shut down pre-maturely, script/download_contig_rerun.py may help. Please read this helper function thoroughly before using it.
#
#
######## Step 4.1: Seperate out Contig ID and WGS ID ######
## To separate contig and wgs
#echo "separating contig and wgs"
#python3 $MY_PATH/source_file/util/extract_contig_wgs_ID.py $WORK_DIR/tmp_contig $contig_list $wgs_list
#wc -l $contig_list $wgs_list
#
#
######## Step 4.2: Download contig file from Contig ID ######
#
## To check on the process finishing situation, helper script: script/last_line_log.py may help.
#ID_list=$contig_list;   outfolder=$CDS_folder/contig_
#num_query=`wc -l $ID_list | cut -d' ' -f1`
#echo "We have $num_query queries in total"
#from=0; to=$num_query;  num_thread=50   # left-close, right-open
#./$MY_PATH/source_file/download_contig.sh $ID_list $outfolder $from $to $num_thread
#wait
#
######## Step 4.3: Download contig file from wgs ID ######
#
#
## For downloading the wgs, the number of thread should be as small as possible. Recommended is 10. Large amount of processes do not work. The process will get shut down by the server.
#ID_list=$wgs_list;              outfolder=$CDS_wgs_folder/wgs_
#num_query=`wc -l $ID_list | cut -d' ' -f1`
#echo "We have $num_query queries in total"
## 16 704
#from=0; to=$num_query;  num_thread=10   # left-close, right-open
#./$MY_PATH/source_file/download_contig.sh $ID_list $outfolder $from $to $num_thread
#wait
#
######## Step 4.4: After processing of all the contig files ######
#rm -r $CDS_folder/temp
#rm -r $CDS_wgs_folder/temp
#mv $CDS_wgs_folder/* $CDS_folder
#
#
######################## Step 5: Extract CDS information from contig file ######
#echo "extracting CDS from raw data"
#cut -f 2,4 $merged_xref_list > $WORK_DIR/tmp_merged_xref_list
#python3 $MY_PATH/source_file/util/extract_CDS_parallel.py $WORK_DIR/tmp_merged_xref_list $CDS_folder $loc_folder 2> $WORK_DIR/extract_CDS_$version.log
#
#
######################## Step 6: Generate id_loc file ######
#echo "generating id_loc file"
#python3 $MY_PATH/source_file/util/genFile.py $loc_folder $id_loc
#python3 $MY_PATH/source_file/util/file_line_redundancy_remove.py $id_loc $TEMP_DIR/id_loc_${version}
#
#
######################## Step 7: Remove all the intermediate files  ######
#
## Remove all the temporary files
#
#rm $WORK_DIR/uniref_id_UNI
#rm $WORK_DIR/uniref_id_UPI
#rm $WORK_DIR/uniref_id_UPI_converted
#rm $WORK_DIR/uniprot_ID_list_redundant
#
#rm -f $WORK_DIR/tmp_contig $WORK_DIR/tmp_coding
#
#rm $WORK_DIR/tmp_merged_xref_list
#
#rm $id_loc

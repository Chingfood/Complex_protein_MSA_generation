#!/bin/bash 

a3m_ffindex=$1
version=$2

#a3m_ffindex=uniclust30_${version}_a3m.ffdata
uniprot_ID_list=uniprot_ID_list_${version}
xref_list_folder=xref_list_${version}
xref_batch_folder=xref_batch_${version}
xref_folder=xref_${version}
merged_xref_list=merged_xref_list_${version}
contig_list=contig_list_${version}
wgs_list=wgs_list_${version}
CDS_folder=CDS_folder_${version}
loc_folder=loc_folder_${version}
id_loc=id_loc_${version}

ID_list=$contig_list;   outfolder=$CDS_folder/contig_
num_query=`wc -l $ID_list | cut -d' ' -f1`
echo "We have $num_query queries in total"
# 9 105 361
from=0; to=$num_query;  num_thread=50   # left-close, right-open
./download_contig.sh $ID_list $outfolder $from $to $num_thread
wait



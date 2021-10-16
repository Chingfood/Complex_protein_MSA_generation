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

num_query=`wc -l $uniprot_ID_list | cut -d' ' -f1`
echo "We have $num_query queries in total"
# 97 884 147 in total
source=UniProtKB/TrEMBL
echo "Downloading $source entries"
#from=086780; to=124173570;  num_thread=10
from=0; to=$num_query;  num_thread=50   # left-close, right-open
#from=0;        to=10000;       num_thread=20   # left-close, right-open
#i=0;   k=7000000;      from=$(($i*$k)); to=$(($(($i+1))*$k))
./download_xref.sh $uniprot_ID_list $source $xref_list_folder $from $to $num_thread



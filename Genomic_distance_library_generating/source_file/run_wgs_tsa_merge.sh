#!/bin/bash

a3m_ffindex=$1
version=$2
uniprot_ID_list=uniprot_ID_list_${version}
xref_list_folder=xref_list_${version}
xref_batch_folder=xref_batch_${version}                                                  xref_folder=xref_${version}
merged_xref_list=merged_xref_list_${version}
contig_list=contig_list_${version}
wgs_list=wgs_list_${version}
CDS_folder=CDS_folder_${version}
CDS_wgs_folder=CDS_wgs_folder_${version}
loc_folder=loc_folder_${version}
id_loc=id_loc_${version}

rm -r $CDS_folder/temp
rm -r $CDS_wgs_folder/temp
mv $CDS_wgs_folder/* $CDS_folder


#!/bin/bash

# ./run.sh uniclust30_2017_10_a3m.ffdata 2017_10
#/mnt/data/conmod_databases/uniclust30_series/uniclust30_2018_08/uniclust30_2018_08_a3m.ffdata
a3m_ffindex=$1
version=$2

#a3m_ffindex=uniclust30_${version}_a3m.ffdata
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

echo "extracting CDS from raw data"
cut -f 2,4 $merged_xref_list > tmp_merged_xref_list
python3 extract_CDS_parallel_new.py tmp_merged_xref_list $CDS_folder $loc_folder 2> extract_CDS_$version.log
rm tmp_merged_xref_list

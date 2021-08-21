#!/bin/bash

# ./run.sh uniclust30_2017_10_a3m.ffdata 2017_10

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

# To verify data class
echo "verifying data"
cut -f 5 $merged_xref_list > tmp_contig
cut -f 4 $merged_xref_list > tmp_coding
wc -l tmp_contig tmp_coding

python3 -c "from enaBrowserTools.python3.utils import is_wgs_sequence, is_sequence
with open('tmp_contig', 'r') as f: assert (is_sequence(_) or is_wgs_sequence(_) for _ in [l[:-1] for l in set(f)])" &

python3 -c "from enaBrowserTools.python3.utils import is_coding
with open('tmp_coding', 'r') as f: assert all(is_coding(_) for _ in f.read().strip().split())" &
wait
echo "verified"

# To separate contig and wgs
echo "separating contig and wgs"
python3 extract_contig_wgs_ID.py tmp_contig $contig_list $wgs_list
wc -l $contig_list $wgs_list
rm -f tmp_contig tmp_coding

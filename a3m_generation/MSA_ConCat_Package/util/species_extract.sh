#!/bin/bash

# ----- usage ------ #
usage()
{
	echo "Version 1.00 [2017-10-01] "
	echo "USAGE: ./species_extract.sh <-i input_fasta> [-c CPU_num] [-d nr] "
	echo "[note]: default CPU_num is 4, nr version is set to 'nr'"
	exit 1
}

if [ $# -lt 1 ];
then
        usage
fi
curdir="$(pwd)"



# ----- get arguments ----- #
#-> optional arguments
cpu_num=4       #-> use 4 CPUs
nr=nr
#-> required arguments 
input_fasta=""

#-> parse arguments
while getopts ":i:c:d:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input_fasta=$OPTARG
		;;
	#-> optional arguments
	c)
		cpu_num=$OPTARG
		;;
	d)
		nr=$OPTARG
		;;
	#-> others
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done


# ------ check required arguments ------ #
if [ ! -f "$input_fasta" ]
then
	echo "input_fasta $input_fasta not found !!" >&2
	exit 1
fi
#-> get job id:
fulnam=`basename $input_fasta`
relnam=${fulnam%.*}


# ---- extract species by NR ---- #
#-> run blast
blast_file=$relnam.blast
ncbi-blast/bin/blastp -num_threads $cpu_num -query $input_fasta -task blastp-fast -db databases/NR/nr -out $blast_file -outfmt "7 qseqid qlen slen qcovhsp pident sseqid staxid bitscore score evalue qstart qend sstart send"
#-> get tax_id
seq_id=`grep -v "#" $blast_file | head -n1 | awk -F"\t" '{print $5}'`
tax_id=`grep -v "#" $blast_file | head -n1 | awk -F"\t" '{print $7}'`
#-> get species_line
species_line=`awk -F"\t" '{print $2}' data/names.dmp_species_taxid_better_fixed3_TaxTree | awk '{if($1==a){print NR}}' a=$tax_id | head -n1`
if [ "$species_line" != "" ]
then
	#-> get species_name
	species_full=`head -n $species_line data/names.dmp_species_taxid_better_fixed3_TaxTree | tail -n1`
	species_name=`head -n $species_line data/names.dmp_species_taxid_better_fixed3_TaxTree | tail -n1 | awk -F"\t" '{print $1}'`
	#-> judge Eukaryotic_or_Prokaryotic
	kingdom_id=`head -n $species_line data/names.dmp_species_taxid_better_fixed3_TaxTree | tail -n1 | awk -F"\t" '{print $NF}' | awk '{print $NF}'`
	kingdom_id_=`echo $kingdom_id | awk -F"(" '{print $1}'`
	kingdom_name=`grep "^$kingdom_id" data/Eukaryotic_or_Prokaryotic | awk -F"\t" '{print $2}' | awk -F"(" '{print $1}' `
	check_name=`grep "^$kingdom_id" data/Eukaryotic_or_Prokaryotic | awk -F"\t" '{print $NF}'`
	#-> echo out
	echo "$check_name|$kingdom_id_|$kingdom_name|$tax_id|$species_name|$seq_id|$species_full"
fi
#-> remove temporary file
#rm -f $blast_file

# ----- success exit ------ #
exit 0


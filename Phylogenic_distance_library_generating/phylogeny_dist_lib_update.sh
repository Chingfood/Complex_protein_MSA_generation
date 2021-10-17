#!/bin/bash

#./run_uniref.sh -t TempDir -r r2d_file -f Uniref.ffdata taxdmp

MY_PATH="`dirname \"$0\"`"

taxdmp=/mnt/home/qingyliu/taxdmp 
r2d=$MY_PATH/Rank_to_Digit
# the directory where it is executed
DIR=`pwd`
# the temp directory used, within $DIR
# omit the -p parameter to create a temporal directory in the default location
TEMP_DIR=`mktemp -d -p "$DIR"`
rm -r $TEMP_DIR

function Usage
{
        echo $0 "[ -t TempDir | -r r2d_file | -f Uniref.ffdata ] taxonomy_database_folder] "
        echo "  This script updates the library for calculating taxonomy distance for MSA. "
        echo "  r2d_file: This is a self defined file that convert the taxonomy rank to digits for after processing. The default location of this file if not specified is located at: taxonomy/Rank_to_Digit"
        echo "  Uniref.ffdata: The sequence database file. Uniref database is used for HHblits. The taxonomy database can be based on all the taxonomy classes of proteins shown up in Uniref library. This is not the database used in complex contact MSA. Normally it is under ffdata extension. An example be: UniRef30_2020_06_a3m.ffdata."
        echo "  TempDir: A temporary directory that store all the intermediate output and final output. If not specified, a random temporary directory will be generated automatically. "
        echo "  taxonomy_database_folder: NCBI taxonomy database downloaded. It contains file like names.dmp and nodes.dmp. The taxonomy library for MSA is built based on all the taxonomy classes shown up in the taxonomy database. "
}


while getopts ":t:r:f:" opt; do
        case ${opt} in
                t )
                  TEMP_DIR=$OPTARG
                  ;;
                r )
                  r2d=$OPTARG
                  ;;
                f )
                  ffdata=$OPTARG
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

if [ $# -ne 1 ]; then
        echo "ERROR: Not all input specified"
        Usage
        exit 1
fi


taxdmp=$1
if [ ! -d $taxdmp ]; then
        echo "ERROR: invalid taxonomy database directory $taxdmp"
        Usage
        exit 1
fi

mkdir -p $TEMP_DIR
echo "Temprory Directory: $TEMP_DIR"
echo "taxdmp directory: $taxdmp"
echo "rd2 file path: $r2d"

if [ ! -z $ffdata ]
then
        echo "ffdata link: $ffdata"
fi

output=$TEMP_DIR

mkdir -p $output/temp

#------ Part I: species in NCBI TaxID -------#

#-> 1.1 dump all species names from names.dmp
awk -F"|" '{print $2}' $taxdmp/names.dmp | awk -F"\t" '{print $2}' > $output/temp/names.dmp_species

python3 $MY_PATH/source_code/original_code/file_line_redundancy_remove.py $output/temp/names.dmp_species $output/temp/names.dmp_species_no
#-> 1.2 dump all species TaxID from names.dmp
awk '{print $1}' $taxdmp/names.dmp > $output/temp/names.dmp_taxid

#-> 1.3 output species and taxid in one file
paste $output/temp/names.dmp_species $output/temp/names.dmp_taxid > $output/temp/names.dmp_species_taxid_better


#-> 2.1 get header line starting with '>' from UniProt20 file
#-> 2.2 extract species in the header line file

if [ ! -z "$ffdata" ]; then
        python3 $MY_PATH/source_code/original_code/OutFile_Species_Uniref.py $ffdata $output/temp/species_out_no
fi

### 3-> 2.3 remove all redundancy of the species from UniProt20
#### This step become redundant because last step has redundancy removed already.
#python3 file_line_redundancy_remove.py $output/temp/species_out $output/temp/species_out_no

#------ Part III: generate taxonomy tree -------#

#-> 3.1 dump all necessary nodes from nodes.dmp
awk '{print $1"\t"$3"\t"$5}' $taxdmp/nodes.dmp > $output/temp/nodes.dmp_taxid_prev_rank

#-> 3.2 to simplify the followed-up analysis, we map the TaxNode (i.e., rank) to digit
#----- file: Rank_to_Digit -----#
#no|1
#superkingdom|2
#kingdom|3
#subkingdom|4
#superphylum|5
#phylum|6
#subphylum|7
#superclass|8
#class|9
#subclass|10
#infraclass|11
#superorder|12
#order|13
#suborder|14
#infraorder|15
#parvorder|16
#superfamily|17
#family|18
#subfamily|19
#tribe|20
#subtribe|21
#genus|22
#subgenus|23
#species|24
#subspecies|25
#varietas|26
#forma|27
#
#-> 3.3 do the mapping of rank to digit
awk '{print $3}' $output/temp/nodes.dmp_taxid_prev_rank > $output/temp/nodes.dmp_rank
cp $output/temp/nodes.dmp_rank $output/temp/nodes.dmp_rank_digit
#--> start grep for each rank
for i in `cat $r2d`; 
do 
	a=`echo $i | cut -d '|' -f 1`; 
	b=`echo $i | cut -d '|' -f 2`; 
	sed "s/^$a$/$b/" $output/temp/nodes.dmp_rank_digit > $output/temp/nodes.dmp_rank_digit_; 
	mv $output/temp/nodes.dmp_rank_digit_ $output/temp/nodes.dmp_rank_digit;
done
rm -f $output/temp/nodes.dmp_rank_digit_ $output/temp/nodes.dmp_rank

#-> 3.4 generate a better form of nodes.dmp_taxid_prev_rank
paste $output/temp/nodes.dmp_taxid_prev_rank $output/temp/nodes.dmp_rank_digit > $output/temp/nodes.dmp_taxid_prev_rank_digit
awk '{print $1"\t"$2"\t"$4"\t"$3}' $output/temp/nodes.dmp_taxid_prev_rank_digit > $output/temp/nodes.dmp_taxid_prev_rank_better
rm -f $output/temp/nodes.dmp_taxid_prev_rank_digit



#==================================================================



#------ Part IV: Trace_Species -------#

#Usage:
#[wangsheng@raptorx6 Taxonomy_ID]$ ./Trace_Species
#Trace_Species <species_taxid_file> <taxid_rank_file> <trace_speices_file>
#
#[note]:
#
##-> 4.1 <species_taxid_file>
#   should be a file containing 'species_string', 'taxid_curnode'
#   e.g., key_input_files/names.dmp_species_taxid_better_fixed3
#   in the following format:
#
#all     1
#root    1
#Bacteria        2
#Monera  2
#Procaryotae     2
#Prokaryota      2
#
#+++++++++++++++++++++

##-> 4.2 <taxid_rank_file>
#   should be a file containing 'taxid_curnode', 'taxid_prevnode', 'rank_digit'
#   e.g., key_input_files/nodes.dmp_taxid_prev_rank_better
#   in the following format:
#
#1       1       1       no
#2       131567  2       superkingdom
#6       335928  22      genus
#7       6       24      species
#9       32199   24      species
#10      1706371 22      genus
#
#+++++++++++++++++++++

#-> 4.3 <trace_speices_file>
#   should be a file only containing 'species_string'


#-> 4.4 running example
$MY_PATH/source_code/executable/Trace_Species $output/temp/names.dmp_species_taxid_better $output/temp/nodes.dmp_taxid_prev_rank_better $output/temp/names.dmp_species_no  1> $output/names.dmp_species_taxid_better_fixed3_TaxTree 2> $output/error_log


if [ ! -z "$ffdata" ]; then
$MY_PATH/source_code/executable/Trace_Species $output/temp/names.dmp_species_taxid_better $output/temp/nodes.dmp_taxid_prev_rank_better $output/temp/species_out_no  1> $output/species_out_no_TaxTree 2> $output/error_log
fi



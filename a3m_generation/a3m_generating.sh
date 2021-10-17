#!/bin/bash

nCPU=20
GPUdevice=0
database=/mnt/local/qingyliu/uniref30_2020_06/UniRef30_2020_06

function Usage
{
        echo $0 "[ -t TempDir ] [-f Uniref_path_header] <seq1> <seq2> [nCPU] [GPU_id]  "
        echo "  This script generate complex MSA for seq1 and seq2."
        echo "  Uniref_path_header: The Uniref library used for HHblits. If not specified, the default Uniref_path_header is: /mnt/local/qingyliu/uniref30_2020_06/UniRef30_2020_06"
        echo "  seq1 and seq2: The sequence of query proteins, usually in fasta format"
        echo "  TempDir: A temporary directory that store all the intermediate output and final output. If not specified, a random temporary directory will be generated automatically. "
	echo "  nCPU: Number of CPU used, default is 20 "
	echo "  GPU_id: The ID of the registered GPU devices used, default is device 0   "
}

# the directory where it is executed
DIR=`pwd`
# the temp directory used, within $DIR
# omit the -p parameter to create a temporal directory in the default location           # remove it first in case user provided a temporary directory
TEMP_DIR=`mktemp -d -p "$DIR"`
rm -r $TEMP_DIR

while getopts ":t:f:" opt; do
        case ${opt} in
                t )
                  TEMP_DIR=$OPTARG
                  ;;
		f )
                  database=$OPTARG
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

if [ $# -lt 2 ]; then
        echo "ERROR: Not all input specified"
        Usage
        exit 1
fi

if [ $# -gt 2 ]; then nCPU=$3; fi
if [ $# -gt 3 ]; then GPUdevice=$4; fi

infile1=`readlink -f $1`
infile2=`readlink -f $2`
#nCPU=$3
#GPUdevice=$4
mkdir -p $TEMP_DIR
echo "Temprory Directory: $TEMP_DIR"

basename1=${infile1%.*}
basename1=${basename1##*/}
basename2=${infile2%.*}
basename2=${basename2##*/}
basename="${basename1}_${basename2}"
basename_="${basename2}_${basename1}"

#source_path=`pwd`/`dirname $0`
source_path="`dirname \"$0\"`"
source_path="`readlink -e $source_path`"
echo "source_path=$source_path"

single_MSA_folder=$TEMP_DIR/Single_MSA
merged_MSA_genomic_folder=$TEMP_DIR/Merged_MSA_Genomic
merged_MSA_phylogenic_folder=$TEMP_DIR/Merged_MSA_Phylogenic
mkdir -p $single_MSA_folder $merged_MSA_genomic_folder $merged_MSA_phylogenic_folder

export HHSUITE="$source_path/MSA_ConCat_Package/hhsuite"
export HHLIB="$source_path/MSA_ConCat_Package/hhsuite"
export MSA_ConCat_Package="$source_path/MSA_ConCat_Package"
export AlignZTM="$source_path/AlignZTM"
export Util_path="$source_path/util/"


peval=20; niter=8
#peval=3; niter=3
eval='1E-'$peval
coverage=-1
id=-1
Meff=.7


#------ export ComplexContactHOME --------$
export ComplexContactHOME="$source_path"




#database=/mnt/data/conmod_databases/uniprot20_series/uniprot20_2016_02/uniprot20_2016_02
#database=/mnt/local/qingyliu/uniref30_2020_06/UniRef30_2020_06
arg4hh="-d $database -e $eval -cpu $nCPU -n $niter -maxfilt 100000000 -diff inf -all -neffmax 20"
$source_path/genMSA.sh $infile1 $single_MSA_folder $coverage "$arg4hh" $Meff &
$source_path/genMSA.sh $infile2 $single_MSA_folder $coverage "$arg4hh" $Meff &
wait



#method='Closest'; coverage=60; identity=90
#method='Closest'; coverage=75; identity=90
method='KM'; coverage=-1; identity=90
min_delta_gene=1
max_delta_gene=20
#GPUdevice=0
nCPU=10
argv="-c $coverage -i $identity -m $min_delta_gene -M $max_delta_gene -j $method -d $GPUdevice -t $nCPU"
$source_path/align_genomic.sh -a $single_MSA_folder/$basename1.a3m -b $single_MSA_folder/$basename2.a3m -o $merged_MSA_genomic_folder/$basename $argv
$source_path/align_genomic.sh -a $single_MSA_folder/$basename2.a3m -b $single_MSA_folder/$basename1.a3m -o $merged_MSA_genomic_folder/$basename_ $argv



#coverage=60; identity=75
coverage=-1; identity=90
#GPUdevice=0
nCPU=10
argv="-c $coverage -i $identity -d $GPUdevice -t $nCPU"
$source_path/align_phylogenic.sh -a $single_MSA_folder/$basename1.a3m_specbloc -b $single_MSA_folder/$basename2.a3m_specbloc -o $merged_MSA_phylogenic_folder/$basename $argv
$source_path/align_phylogenic.sh -a $single_MSA_folder/$basename2.a3m_specbloc -b $single_MSA_folder/$basename1.a3m_specbloc -o $merged_MSA_phylogenic_folder/$basename_ $argv



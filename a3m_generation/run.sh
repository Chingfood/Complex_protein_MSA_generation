#!/bin/bash

nCPU=20
GPUdevice=0

if [ $# -lt 2 ]; then echo "Usage: ./run.sh <seq1> <seq2> [nCPU] [GPU_id]"; exit; fi
if [ $# -gt 2 ]; then nCPU=$3; fi
if [ $# -gt 3 ]; then GPUdevice=$4; fi

infile1=`readlink -f $1`
infile2=`readlink -f $2`
#nCPU=$3
#GPUdevice=$4

basename1=${infile1%.*}
basename1=${basename1##*/}
basename2=${infile2%.*}
basename2=${basename2##*/}
basename="${basename1}_${basename2}"
basename_="${basename2}_${basename1}"

source_path=`pwd`/`dirname $0`
echo "source_path=$source_path"

single_MSA_folder="Single_MSA"
merged_MSA_genomic_folder="Merged_MSA_Genomic"
merged_MSA_phylogenic_folder="Merged_MSA_Phylogenic"
prediction_genomic_folder="Prediction_Genomic"
prediction_phylogenic_folder="Prediction_Phylogenic"
prediction_merged_folder="Prediction_Merged"
mkdir -p $single_MSA_folder $merged_MSA_genomic_folder $merged_MSA_phylogenic_folder $prediction_genomic_folder $prediction_phylogenic_folder $prediction_merged_folder

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
database=/mnt/local/qingyliu/uniref30_2020_06/UniRef30_2020_06
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



#GPUdevice=0
function f {
#	for extension in seq a3m a2m tgt ss3 acc pot ccmpred; do cp $1/$4.$extension.tar.gz $2/; done
	for extension in seq a3m a2m tgt ss3 acc pot; do cp $1/$4.$extension.tar.gz $2/; done
	$source_path/predict.sh $1/$4.a3m $2 $3 $GPUdevice
}
f $merged_MSA_genomic_folder $prediction_genomic_folder genomic $basename
f $merged_MSA_phylogenic_folder $prediction_phylogenic_folder phylogenic $basename
f $merged_MSA_genomic_folder $prediction_genomic_folder genomic $basename_
f $merged_MSA_phylogenic_folder $prediction_phylogenic_folder phylogenic $basename_
wait



#---- merge predicted interfacial contact map -----#
for folder in $prediction_genomic_folder $prediction_phylogenic_folder; do python $source_path/util/swap.py $infile1 $infile2 $folder/$basename_.gcnn    $folder/$basename_.gcnn.swap; done
#for folder in $merged_MSA_genomic_folder $merged_MSA_phylogenic_folder; do python $source_path/util/swap.py $infile1 $infile2 $folder/$basename_.ccmpred $folder/$basename_.ccmpred.swap; done
python $source_path/util/merge.py $infile1 $infile2 $prediction_genomic_folder/$basename.gcnn    $prediction_phylogenic_folder/$basename.gcnn    $prediction_genomic_folder/$basename_.gcnn.swap    $prediction_phylogenic_folder/$basename_.gcnn.swap    $prediction_merged_folder/$basename.gcnn_
#python $source_path/util/merge.py $infile1 $infile2 $merged_MSA_genomic_folder/$basename.ccmpred $merged_MSA_phylogenic_folder/$basename.ccmpred $merged_MSA_genomic_folder/$basename_.ccmpred.swap $merged_MSA_phylogenic_folder/$basename_.ccmpred.swap $prediction_merged_folder/$basename.ccmpred_

#-> swap and merge epad_prob predictions
for folder in $prediction_genomic_folder $prediction_phylogenic_folder; do python $source_path/util/swapEPAD.py $infile1 $infile2 $prediction_genomic_folder/$basename_.epad_prob $folder/$basename_.epad_prop.swap; done
python $source_path/util/mergeEPAD.py $infile1 $infile2 $prediction_genomic_folder/$basename.epad_prob $prediction_phylogenic_folder/$basename.epad_prob $prediction_genomic_folder/$basename_.epad_prop.swap $prediction_phylogenic_folder/$basename_.epad_prop.swap $prediction_merged_folder/$basename.epad_



#---- draw a heatmap of the predicted interfacial contact map -----#
len1=`tail -n1 $single_MSA_folder/$basename1.seq | wc | awk '{print $3-1}'`
len2=`tail -n1 $single_MSA_folder/$basename2.seq | wc | awk '{print $3-1}'`
python $source_path/util/list2heatmap.py $prediction_merged_folder/$basename.gcnn_inter_list $len1 $len2 $prediction_merged_folder/$basename.gcnn_inter.png
#python $source_path/util/list2heatmap.py $prediction_merged_folder/$basename.ccmpred_inter_list $len1 $len2 $prediction_merged_folder/$basename.ccmpred_inter.png


#---- generate top K list ----- #
alpha=5
topk=$(( alpha*(len1+len2) ))
head -n $topk $prediction_merged_folder/$basename.gcnn_inter_list | awk '{print $0}' > $prediction_merged_folder/$basename.gcnn_inter_list_top
#head -n $topk $prediction_merged_folder/$basename.ccmpred_inter_list | awk '{print $0}' > $prediction_merged_folder/$basename.ccmpred_inter_list_top

#-> grep topk epad_prob predictions
awk '{print $1" "$2" "}' $prediction_merged_folder/$basename.gcnn_inter_list_top > $prediction_merged_folder/$basename.gcnn_inter_list_top_
rm -f $prediction_merged_folder/$basename.epad_inter_top
len=`wc $prediction_merged_folder/$basename.gcnn_inter_list_top_ | awk '{print $1}'`
for ((i=1;i<=$len;i++))
do
	content=`head -n $i $prediction_merged_folder/$basename.gcnn_inter_list_top_ | tail -n1`
	grep "^$content" $prediction_merged_folder/$basename.epad_inter >> $prediction_merged_folder/$basename.epad_inter_top
done
rm -f $prediction_merged_folder/$basename.gcnn_inter_list_top_


#---- put all required files to a zip folder -------#
#-> simplified zip for user
mkdir -p ${basename}.for_user
#-> copy readme
cp $Util_path/0README_simp ${basename}.for_user
#-> original sequence in FASTA format
cp $infile1 ${basename}.for_user/${basename1}.fasta
cp $infile2 ${basename}.for_user/${basename2}.fasta
#-> final merged contact prediction
if [  ! -s "$prediction_merged_folder/$basename.gcnn_inter" ]
then
	echo "FAILED !! $prediction_merged_folder/$basename.gcnn_inter not generated."
	exit 1
else
	cp $prediction_merged_folder/$basename.gcnn_inter ${basename}.for_user
fi
if [ ! -s "$prediction_merged_folder/$basename.gcnn_inter.png" ]
then
	echo "FAILED !! $prediction_merged_folder/$basename.gcnn_inter.png not generated."
	exit 1
else
	cp $prediction_merged_folder/$basename.gcnn_inter.png ${basename}.for_user
fi
if [ ! -s "$prediction_merged_folder/$basename.gcnn_inter_list_top" ]
then
	echo "FAILED !! $prediction_merged_folder/$basename.gcnn_inter_list_top not generated."
	exit 1
else
	cp $prediction_merged_folder/$basename.gcnn_inter_list_top ${basename}.for_user/$basename.gcnn_inter_top
fi
if [ ! -s "$prediction_merged_folder/$basename.epad_inter" ]
then
	echo "FAILED !! $prediction_merged_folder/$basename.epad_inter not generated."
	exit 1
else
	cp $prediction_merged_folder/$basename.epad_inter ${basename}.for_user
fi
if [ ! -s "$prediction_merged_folder/$basename.epad_inter_top" ]
then
	echo "FAILED !! $prediction_merged_folder/$basename.epad_inter_top not generated."
	exit 1
else
	cp $prediction_merged_folder/$basename.epad_inter_top ${basename}.for_user
fi
#-> zip 
zip -r ${basename}.for_user.zip ${basename}.for_user/


#----------- detailed file for download ---------------#
zip=1    #-> 0 for simp version; 1 for full version
mkdir -p ${basename}.all_in_one
#-> original sequence in FASTA format
cp $infile1 ${basename}.all_in_one/${basename1}.fasta
cp $infile2 ${basename}.all_in_one/${basename2}.fasta
#-> final merged contact prediction
cp $prediction_merged_folder/$basename.gcnn_inter* ${basename}.all_in_one
cp $prediction_merged_folder/$basename.epad_inter* ${basename}.all_in_one
#-> copy others
if [ $zip -eq 0 ]
then
	#--| copy readme
	cp $Util_path/0README_simp ${basename}.all_in_one
	rm -rf ${basename}.all_in_one/$basename.gcnn_inter_list
	mv ${basename}.all_in_one/$basename.gcnn_inter_list_top ${basename}.all_in_one/$basename.gcnn_inter_top
else
	#--| copy readme
	cp $Util_path/0README_total ${basename}.all_in_one
	#-> single MSA in A3M format
	cp $single_MSA_folder/${basename1}.a3m ${basename}.all_in_one
	cp $single_MSA_folder/${basename2}.a3m ${basename}.all_in_one
	#-> paired MSA in A3M format
	cp $merged_MSA_genomic_folder/${basename}.a3m ${basename}.all_in_one/${basename}.a3m_genomic
	cp $merged_MSA_phylogenic_folder/${basename}.a3m ${basename}.all_in_one/${basename}.a3m_phylogenic
	#-> contact prediction by single method
	cp $prediction_genomic_folder/${basename}.gcnn ${basename}.all_in_one/${basename}.gcnn_genomic
	cp $prediction_phylogenic_folder/${basename}.gcnn ${basename}.all_in_one/${basename}.gcnn_phylogenic
	cp $prediction_genomic_folder/${basename}.epad_prob ${basename}.all_in_one/${basename}.epad_prob_genomic
	cp $prediction_phylogenic_folder/${basename}.epad_prob ${basename}.all_in_one/${basename}.epad_prob_phylogenic
	#-> remove top
	rm -f ${basename}.all_in_one/$basename.gcnn_inter_list_top
	#-> copy intra predictions
	cp $prediction_merged_folder/$basename.gcnn_intra_? ${basename}.all_in_one
	cp $prediction_merged_folder/$basename.epad_intra_? ${basename}.all_in_one
fi

#--- make a zip package ---#
zip -r ${basename}.all_in_one.zip ${basename}.all_in_one/

#---------- exit --------------#
exit 0


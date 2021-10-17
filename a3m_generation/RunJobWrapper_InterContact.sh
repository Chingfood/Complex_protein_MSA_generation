#!/bin/bash
#$ -S /bin/bash

#------- main path --------#
CMD_DIR=~/GitBucket/InterContact_Pred
CPU_NUMBER=20
GPU_ID=0
#-> remote machine
#LOCALE_MAC=""             #-> set "" for using default local machine, such as '$USER@hostname`
#REMOTE_MAC=""             #-> set "" for not using remote machine to run GCNN for contact prediction
#LOCALE_MAC=RaptorX@raptorx3.uchicago.edu   #-> set raptorx3 as local machine to run basic InterContact_Pred
#REMOTE_MAC=RaptorX@raptorx6.uchicago.edu   #-> set raptorx6 as remote machine to run GCNN for contact prediction
export LOCALE_MAC
export REMOTE_MAC

#------- usage ------------#
if [ $# -lt 3 ]
then
        echo "Usage: ./RunJobWrapper_InterContact.sh resultLocation seq_file1 seq_file2 [email] [qsub] "
	echo "       resultLocation shall be something like RaptorX@raptorx:resultDir or RaptorX@localhost:resultDir or localhost:resultDir "
	echo "       seq_file looks like XXXX.seq or XXX.fasta and is saved at resultLocation "
	echo "       [qsub] mode will only be applied on cruncher machine, with job names as 'CrunQsub' "
        exit 1
fi

#------ argument ----------#
resultLocation=$1
seqFile1=$2
seqFile2=$3
email=""
if [ $# -gt 4 ]
then
	email=$5
fi

#------ addi command for qsub -----#
if [ $# -gt 5 ]
then
	CPU_NUMBER=$NSLOTS
fi

#----- sequence name -------# 
seqNumber1=`echo $seqFile1 | cut -f1 -d'.'`
fulnam1=`basename $seqFile1`
relnam1=${fulnam1%.*}
seqNumber2=`echo $seqFile2 | cut -f1 -d'.'`
fulnam2=`basename $seqFile2`
relnam2=${fulnam2%.*}

#-> tmp folder
tmpdir="${relnam1}_${relnam2}"


#---------------- running job --------------#
## Change directory to the work direcotry for the property prediction server 
cd ${CMD_DIR}
mkdir -p $tmpdir
cd $tmpdir

## obtain the sequence file by scp
scp ${resultLocation}/${seqFile1} ./$relnam1.fasta
rcode=$?
if [ $rcode -ne 0 ]
then
	echo "Failed to scp the input sequence file 1 from the remote folder ${resultLocation}"
	exit 1
fi
scp ${resultLocation}/${seqFile2} ./$relnam2.fasta
rcode=$?
if [ $rcode -ne 0 ]
then
	echo "Failed to scp the input sequence file 2 from the remote folder ${resultLocation}"
	exit 1
fi

## running a property prediction job
echo "running a interfacial contact prediciton job on ${HOSTNAME} for $seqFile1 and $seqFile2 ..."
real_cmd="../run.sh"
echo "${real_cmd} $relnam1.fasta $relnam2.fasta ${CPU_NUMBER} ${GPU_ID} "
${real_cmd} $relnam1.fasta $relnam2.fasta ${CPU_NUMBER} ${GPU_ID}
rcode=$?
if [ $rcode -ne 0 ]
then
	echo "Failed to execute ${real_cmd} $relnam1.fasta $relnam2.fasta ${CPU_NUMBER} ${GPU_ID} "
	exit 1
fi

## remove the sequence file by scp 
cd ../

#---------------- copy/remove --------------#
##copy the results to the resultLocation
scp -r ${CMD_DIR}/${tmpdir}/* ${resultLocation}/
rcode=$?
if [ $rcode -ne 0 ]
then
	echo "Failed to scp the prediction results back to ${resultLocation}"
	exit 1
fi

## remove temporary files
rm -r ${CMD_DIR}/${tmpdir}


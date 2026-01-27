#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=20:mem=500gb:scratch_local=1000gb:spec=7.0
#PBS -j oe
#PBS -N dDocent_fullRun
#PBS -m abe
#--------------------------------------------------------------------------------------------

#Run dDocent with all samples (after optimizing parameters - similarity, k1, k2)

server=brno12-cerit
cores=$TORQUE_RESC_TOTAL_PROCS
#define dir with fq.gz files
DATADIR="/storage/brno12-cerit/home/tomasfer/louky/TACR/Galium2025/fastq/demultiplexed/subsample1mil"
workdir=results
#workdir=$(sed "s/.*\///" <<< ${DATADIR})
#move to SCRATCHDIR
cd $SCRATCHDIR

#copy data from datadir
mkdir ${workdir}
cp ${DATADIR}/*.gz .

#run cutadapt first (remove bases of restriction sites) and save the modified files to 'workdir'
#activate cutadapt
mamba activate cutadapt-5.2
#this works for Walter's coding scheme
#indnames=$(find . -type f | grep -Eo '[A-Z]{3}\.[A-Z]{3}.*[A-Z][0-9]{2}' | sort | uniq)
#this works for our samples
indnames=$(ls *.fq.gz | cut -d'.' -f1)
for indname in $indnames; do
	fastq1="./${indname}.R.fq.gz"
	fastq2="./${indname}.F.fq.gz"
	# --cores=0 auto-detects number of available cores
	cutadapt --cores=0 --cut 3 -U 5 --output ${workdir}/${fastq1} --paired-output ${workdir}/${fastq2} ${fastq1} ${fastq2}
done
mamba deactivate

#run dDocent
cd ${workdir}
cp ${DATADIR}/config.file .
ls *.gz > inputgzlist
#create dDocent instalation with mamba (just once, then comment)
#mamba create --prefix /storage/brno2/home/${LOGNAME}/test_env_ddoc ddocent
mamba activate /auto/brno2/home/${LOGNAME}/test_env_ddoc
dDocent config.file
mamba deactivate
#delete input data in workdir
for i in $(cat inputgzlist); do
	rm ${i}
done

#copy data back home
cd ..
cp -r ${workdir} ${DATADIR}

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi

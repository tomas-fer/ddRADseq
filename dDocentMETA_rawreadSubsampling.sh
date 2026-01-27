#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=256gb:scratch_local=200gb:spec=7.0
#PBS -j oe
#PBS -N dDocent_rawreadsSubsampling
#PBS -m abe
#--------------------------------------------------------------------------------------------

#raw read subsampling and renaming

server=brno12-cerit
cores=$TORQUE_RESC_TOTAL_PROCS
#define dir with demultiplexed files (*.{1,2}.fq.gz)
DATADIR="/storage/brno12-cerit/home/${LOGNAME}/louky/TACR/Galium2025/fastq/demultiplexed"
#subsample to nr. read pairs
subsam=1000000

#subsample using seqtk
echo -e "Subsampling reads to $subsam..."
cd ${DATADIR}
mkdir $subsam
for file in $(ls *.fq.gz | cut -d'.' -f1,2); do
	echo ${file}
	seqtk sample -s100 ${file}.fq.gz ${subsam} > ${subsam}/${file}.fq
	gzip ${subsam}/${file}.fq
done

#rename samples to follow dDocent naming scheme (from .{1,2}.fq.gz to .{F,R}.fq.gz
echo -e "Renaming samples to {F,R}..."
cd $subsam
indnames=$(ls *.fq.gz | cut -d'.' -f1)
for indname in $indnames; do
	mv ${indname}.1.fq.gz ${indname}.F.fq.gz
	mv ${indname}.2.fq.gz ${indname}.R.fq.gz
done
echo -e "Finished..."

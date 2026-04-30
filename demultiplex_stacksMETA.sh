#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=2:ncpus=4:mem=400gb:scratch_local=500gb
#PBS -j oe
#PBS -N demultiplex_stacks
#PBS -m abe

#data
folder=/storage/brno12-cerit/home/tomasfer/louky/TACR/Galium2025
folder=${folder}/fastq
file1=GIUS-AAMM-0001-99AD-1N1_1.fastq.gz
file2=GIUS-AAMM-0001-99AD-1N1_2.fastq.gz
barcodes=barcodes_Galium1.txt

#MetaCentrum modules
module add stacks/2.68-gcc-10.2.1-hqulvik
module add jdk-8
module add seqtk-1.0
module add vcftools-0.1.16
module add python36-modules-gcc

#start
cd $SCRATCHDIR

#copy files
echo "Copying input files..."
cp ${folder}/${file1} .
cp ${folder}/${file2} .
cp ${folder}/${barcodes} .

#demultiplex original files
echo "Demultiplexing (process_radtags)..."
mkdir demultiplexed
#process_radtags -1 ${file1} -2 ${file2} -b ${barcodes} -o demultiplexed --renz_1 ecoRI --renz_2 mspI -c -q
process_radtags --paired -1 ${file1} -2 ${file2} -o demultiplexed -i gzfastq -y gzfastq -b ${barcodes} --barcode_dist_1 2 --barcode_dist_2 2 --inline_index --renz_1 'ecoRI' --renz_2 'mspI' --disable_rad_check --retain_header


echo "Copying demultiplexed data back..."
cp -r demultiplexed ${folder}

echo "Cuting out first 3 or 5 bases (cutadapt)..."
module add py-cutadapt
cd demultiplexed
mkdir cutadapt
indnames=$(ls *1.fq.gz | sed 's/\.1\.fq\.gz//')

#CORRECT version is cutadapt --cut 3 -U 5 !!!!
for i in $indnames; do
	fastq1="./${i}.R.fq.gz"
	fastq2="./${i}.F.fq.gz"
	cutadapt --cores=$TORQUE_RESC_TOTAL_PROCS --cut 3 -U 5 --output "cutadapt/${fastq1}" --paired-output "cutadapt/${fastq2}" $fastq1 $fastq2
done

echo "Copying cutadapt data back..."
cp -r cutadapt ${folder}

echo "Done..."
exit

#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=2:ncpus=4:mem=400gb:scratch_local=500gb
#PBS -j oe
#PBS -N demultiplex_stacks
#PBS -m abe

#data
folder=/storage/brno12-cerit/home/tomasfer/louky/TACR/Galium2025/fastq
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
process_radtags -1 ${file1} -2 ${file2} -b ${barcodes} -o demultiplexed --renz_1 ecoRI --renz_2 mspI -c -q

echo "Copying demultiplexed data back..."
cp -r demultiplexed ${folder}

echo "Done..."
exit

#trim data with Trimmomatic (10bp at the beginning of each fq.gz)
mkdir trimmed

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

for file in $(ls demultiplexed_GCCAAT/*.fq.gz | grep -v rem | cut -d'/' -f2 | cut -d'.' -f1,2); do echo ${file}; java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 6 -phred33 demultiplexed_GCCAAT/${file}.fq.gz trimmed/${file}.trimmed.fq.gz HEADCROP:10; done

for file in $(ls demultiplexed_CTTGTA/*.fq.gz | grep -v rem | cut -d'/' -f2 | cut -d'.' -f1,2); do echo ${file}; java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 6 -phred33 demultiplexed_GCCAAT/${file}.fq.gz trimmed/${file}.trimmed.fq.gz HEADCROP:10; done

#subsample data to 1mil. reads
mkdir trimmed/subsample1mil

for file in $(ls trimmed/*.fq.gz | cut -d'/' -f2 | cut -d'.' -f1,2); do echo ${file}; seqtk sample -s100 trimmed/${file}.trimmed.fq.gz 1000000 > trimmed/subsample1mil/${file}.trimmed.sub1mil.fq; gzip trimmed/subsample1mil/${file}.trimmed.sub1mil.fq; done

#sort data according species
species=Loricaria
mkdir trimmed/subsample1mil/${species}
mv trimmed/subsample1mil/LR* trimmed/subsample1mil/${species}

#rename trimmed samples
for file in $(ls trimmed/subsample1mil/${species}/*.fq.gz | cut -d'/' -f4 | cut -d'.' -f1,2); do mv ${file}.trimmed.sub1mil.fq.gz ${file}.fq.gz; done

#prepare list of samples
ls trimmed/subsample1mil/${species}/*.fq.gz | cut -d'/' -f4 | cut -d'.' -f1 | sort | uniq > samples.txt
#use this samples.txt for preparation of popmap file and save it to trimmed/subsample1mil/${species}/

#run Stacks de novo map
mkdir -p STACKS_denovo_map/batch_1_m5_n4_M4
inputdir=trimmed/subsample1mil/${species}
popmap=trimmed/subsample1mil/${species}/popmap_${species}.txt
outputdir=STACKS_denovo_map/batch_1_m5_n4_M4

denovo_map.pl -m 5 -n 4 -M 4 -X "ustacks:-H" --samples $inputdir --popmap $popmap -o $outputdir -T 6 --paired

#optimization of '-M' and '-n' (M=n in our case)
#loop with 8 individuals


#populations
inputdir=STACKS_denovo_map/batch_1_m5_n4_M4

#RAxML output
#popmap - everything is single population for phylogenetics
#--max-obs-het 0.65 (only loci with maximum heterozygosity of 0.65)
outfiles=STACKS_denovo_map/batch_1_m5_n4_M4/exportRAxML
popmap=trimmed/subsample1mil/${species}/popmap_${species}_1pop.txt
populations -P $inputdir -O $outfiles -M $popmap -t 4 --max-obs-het 0.65 --vcf
#filtering for missing data etc.
vcf=STACKS_denovo_map/batch_1_m5_n4_M4/exportRAxML/populations.snps.vcf
#Keep only SNPs present in at least 20%‚ 50%, 70% indivs
vcftools --vcf $vcf --remove-indels --max-missing 0.2 --minDP 10 --out ${outfiles}/${species}.max20percentMissing_minDP10 --recode
#After filtering, kept 133036 out of a possible 293786 Sites
vcftools --vcf $vcf --remove-indels --max-missing 0.5 --minDP 10 --out ${outfiles}/${species}.max50percentMissing_minDP10 --recode
#After filtering, kept 28029 out of a possible 293786 Sites
vcftools --vcf $vcf --remove-indels --max-missing 0.7 --minDP 10 --out ${outfiles}/${species}.max70percentMissing_minDP10 --recode
#After filtering, kept 9089 out of a possible 293786 Sites

wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
python vcf2phylip.py -i ${outfiles}/$species.max70percentMissing_minDP10.recode.vcf -n -f -b

#-p 2 (SNP in at least 2 pops)
#-r 0.5 (SNP in at least 50% of indivs)
#--write_single_snp (only write 1st SNP per locus)
#popmap - 3 populations
outfiles=STACKS_denovo_map/batch_1_m5_n4_M4/exportSTRUCTURE
popmap=trimmed/subsample1mil/Loricaria/popmap_Loricaria.txt
populations -P $inputdir -O $outfiles -M $popmap -p 2 -r 0.5 -t 4 --max-obs-het 0.65 --write-single-snp --vcf --structure

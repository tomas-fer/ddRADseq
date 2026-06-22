#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=2gb
#PBS -j oe
#PBS -N RpackagesSambaR_setup
#PBS -m abe

#define variables
server=brno12-cerit
datafolder=louky/TACR/Lychnis2026/fastq/subsamp1000000/results/filtered.vcf
species=LYC.FLO
name=${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode

echo "This runs LEA in SambaR for:"
echo "Species: $species"
echo "File: ${name}.vcf"
echo -e "Folder: /storage/${server}/home/${LOGNAME}/${datafolder}\n"

cd $SCRATCHDIR

#set modules
module add gdal/3.4.3-gcc-10.2.1-p5yz4uq
module add proj/8.2.1-gcc-10.2.1-iwko2my
module add udunits/2.2.28-gcc-10.2.1-36jrazx
module add vcftools
module add r/4.4.0-gcc-10.2.1-ssuwpvb

export R_LIBS="/storage/${server}/home/${LOGNAME}/RpackagesSambaR"

#download plink1.9
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
unzip plink_linux_x86_64_20250819.zip

#copy SambaR files
wget https://raw.githubusercontent.com/mennodejong1986/SambaR/refs/heads/master/SAMBAR_v1.10.txt
#cp /storage/${server}/home/${LOGNAME}/RpackagesSambaR/SAMBAR_v1.10.txt .
cp /storage/${server}/home/${LOGNAME}/RpackagesSambaR/mypackageslist.txt .
cp /storage/${server}/home/${LOGNAME}/RpackagesSambaR/runLEAinSambaR.R .

#copy data
cp /storage/${server}/home/${LOGNAME}/${datafolder}/${name}.vcf .

#convert data
echo -e "\nConverting data...\n"
gzip ${name}.vcf
vcftools --gzvcf ${name}.vcf.gz --plink --out ${name}
#correct PED file - 1st column should contain pop name (taken from the full name 'popname_samplename' in the second column)
cut -f1 ${name}.ped | cut -f1 -d'_' > pops
cut -f2- ${name}.ped > rest
paste pops rest > ${name}.ped
rm pops rest
#calculate depth
echo -e "\nCalculating depth...\n"
vcftools --gzvcf ${name}.vcf.gz --depth
vcftools --gzvcf ${name}.vcf.gz --site-mean-depth
#Convert PED/MAP to RAW/BIM
./plink --file ${name} --chr-set 95 --allow-extra-chr --make-bed --recode A --out ${name}

#Run LEA in SambaR
echo -e "\nRunning SambaR...\n"
R --slave -f runLEAinSambaR.R $name

#Plot maps from LEA results
echo -e "\nPlotting maps...\n"
cd SambaR_output/Structure
cp /storage/${server}/home/${LOGNAME}/${datafolder}/localitiesCoor.txt .
wget https://raw.githubusercontent.com/tomas-fer/ddRADseq/refs/heads/main/LEAmakePieMaps.sh
wget https://raw.githubusercontent.com/tomas-fer/ddRADseq/refs/heads/main/plotLEA_maps.R
chmod +x LEAmakePieMaps.sh
chmod +x plotLEA_maps.R
./LEAmakePieMaps.sh
cd ../..

#Copy results home
mkdir /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
cp out.* /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
cp pairwise_missingness.txt /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
rm -r SambaR_output/{Demography,Divergence,Diversity,Inputfiles,Kinship,Maps,Selection}
cp -r SambaR_output /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults

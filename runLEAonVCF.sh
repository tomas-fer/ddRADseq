#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=4gb:scratch_local=2gb
#PBS -j oe
#PBS -N runLEAonVCF
#PBS -m abe

#define variables
server=brno12-cerit
datafolder=louky/TACR/Lychnis2026/fastq/subsamp1000000/results/filtered.vcf
species=LYC.FLO
fullname="Lychnis flos-cuculi"
name=${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode
K=10

echo "This runs LEA in SambaR for:"
echo "Species full name: $fullname"
echo "Species code: $species"
echo "File: ${name}.vcf"
echo "K: ${K}"
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
R --slave -f runLEAinSambaR.R $name $K

#Plot maps from LEA results
echo -e "\nPlotting maps...\n"
cd SambaR_output/Structure
cp /storage/${server}/home/${LOGNAME}/${datafolder}/localitiesCoor.txt .
wget https://raw.githubusercontent.com/tomas-fer/ddRADseq/refs/heads/main/LEAmakePieMaps.sh
wget https://raw.githubusercontent.com/tomas-fer/ddRADseq/refs/heads/main/plotLEA_maps.R
chmod +x LEAmakePieMaps.sh
chmod +x plotLEA_maps.R
./LEAmakePieMaps.sh $K $species

#Plot ancestry visualization
echo -e "\nPlotting ancestry maps...\n"
cp /storage/${server}/home/${LOGNAME}/RpackagesSambaR/POPSutilitiesCorrectedConstraines.r .
cp /storage/${server}/home/${LOGNAME}/RpackagesSambaR/Europe.asc .
#Modify file with locality coordinates
cut -f2,3 localitiesCoor.txt | tail -n +2 > coordinates.txt
for i in $(seq 2 $K); do
	cut -f2- Structureplot.LEAqmatrix_K${i}.txt > k${i}.Q
done
echo ${fullname} > fullname
R --slave -f PlotMapFromQmatrix.R $K
#Merge all PDFs into a single file
echo -e"\nMerging PDFs\n"
#Check if there is pdfbox.jar available
if [ -f "./pdfbox.jar" ]; then
	echo "PDFbox found..."
else
	echo "Downloading PDFbox..."
	#Check the newest PDFbox version and silently download it
	pdfboxver=$(wget -q -O- https://downloads.apache.org/pdfbox/ | grep "2\.0\." | cut -d'"' -f6 | sed 's/.$//')
	#download the newest version and rename
	wget -q https://downloads.apache.org/pdfbox/${pdfboxver}/pdfbox-app-${pdfboxver}.jar
	mv pdfbox-app-${pdfboxver}.jar pdfbox.jar
fi
#merge K*.pdf to a single PDF
java -jar pdfbox.jar PDFMerger K*_regions.pdf ${species}_K2-${k}_regions.pdf

cd ../..

#Copy results home
mkdir /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
cp out.* /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
cp pairwise_missingness.txt /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults
rm -r SambaR_output/{Demography,Divergence,Diversity,Inputfiles,Kinship,Maps,Selection}
cp -r SambaR_output /storage/${server}/home/${LOGNAME}/${datafolder}/LEAresults

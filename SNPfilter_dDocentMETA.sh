#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=256gb:scratch_local=200gb:spec=7.0
#PBS -j oe
#PBS -N dDocent_SNPfilter
#PBS -m abe
#--------------------------------------------------------------------------------------------

#SNP filtering after dDocent full run

server=brno12-cerit
cores=$TORQUE_RESC_TOTAL_PROCS
#define dir with VCF input files (TotalRawSNPs.vcf)
DATADIR="/storage/brno12-cerit/home/${LOGNAME}/louky/TACR/Galium2025/fastq/demultiplexed/subsample1mil/results"
species=GAL.MOL

mkdir ${DATADIR}/filtered.vcf
#move to SCRATCHDIR
cd $SCRATCHDIR
mkdir ${species}
cd ${species}

#copy input VCF
echo -e "\nCopying data..."
cp ${DATADIR}/TotalRawSNPs.vcf .

#activate dDocent
echo -e "\nActivating dDocent..."
mamba activate /auto/brno2/home/${LOGNAME}/test_env_ddoc

#decompose complex variants into SNPs and others
echo -e "\nDecomposing complex variants (vcfallelicprimitives)..."
vcfallelicprimitives -kg TotalRawSNPs.vcf > "$species".TotalRawSNPsPrim.vcf

#keep variants (SNPs), i.e. remove INDELs
#remove loci not genotyped in >50% of individuals, min quality score=30 minor allele count of 3
echo -e "\nRemoving loci not genotyped in >50% of indiv..."
vcftools --gzvcf "$species".TotalRawSNPsPrim.vcf --remove-indels --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out "$species".SNPs.g5mac3
#keep genotypes with >= 3 reads and only biallelic SNPs
echo -e "\nKeep genotypes with >= 3 reads and only biallelic SNPs..."
vcftools --vcf "$species".SNPs.g5mac3.recode.vcf --minDP 3 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out "$species".SNPs.g5mac3minDP3bi
#output observed and expected heterozygosities and statistics
echo -e "\nOutput Ho and He heterozygosities..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3bi.recode.vcf --hardy
#search for loci that have an excess of heterozygotes P_HET_EXCESS <1e-15 - to exclude loci that are fixed heterozygotes
echo -e "\nExclude fixed heterozygous loci..."
cat 'out.hwe' | awk '{if($8<1e-15){print $0}}' | cut -f1-2 > het_excess15.loci
#remove them
vcftools --gzvcf "$species".SNPs.g5mac3minDP3bi.recode.vcf --exclude-positions het_excess15.loci --recode --recode-INFO-all --out "$species".SNPs.g5mac3minDP3birl15
#output a list of all individuals with the percentage of missing genotypes (default filename "out.imiss") - and remove indiv with too missing values
echo -e "\nRemove indiv with too many missing data..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15.recode.vcf --missing-indv
#search for individuals that have > X missing genotypes
cat 'out.imiss' | awk '{if($5>0.75){print $0}}' | cut -f1 > lowDP75.indv
#remove them
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15.recode.vcf --remove lowDP75.indv --recode --recode-INFO-all --out "$species".SNPs.g5mac3minDP3birl15rm75
#filter loci by allele balance - these are from dDocent manual - check!
echo -e "\nFilter loci by allele balance (F1)..."
vcffilter -s -f  "AB > 0.2 & AB < 0.8 | AB < 0.01 | AB > 0.99" "$species".SNPs.g5mac3minDP3birl15rm75.recode.vcf > "$species".SNPs.g5mac3minDP3birl15rm75.recode.F1.vcf
#filter out SNPs that have reads from both stands
echo -e "\nFilter out SNPs that have reads from both stands (F2)..."
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s "$species".SNPs.g5mac3minDP3birl15rm75.recode.F1.vcf > "$species".SNPs.g5mac3minDP3birl15rm75.recode.F2.vcf
#filter out loci hat have strongly different mapping qualities of two alleles
echo -e "\nFilter out loci hat have strongly different mapping qualities of two alleles (F3)..."
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" "$species".SNPs.g5mac3minDP3birl15rm75.recode.F2.vcf > "$species".SNPs.g5mac3minDP3birl15rm75.recode.F3.vcf
#filter out loci if ref and alt allels are differently paired (ref = paired; alt = non paired)
echo -e "\nFilter out loci if ref and alt allels are differently paired (F4)..."
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s "$species".SNPs.g5mac3minDP3birl15rm75.recode.F3.vcf > "$species".SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf
#keep SNPs genotyped in >66% of individuals; mean read depth for genotype call = 5/10/20
echo -e "\nKeep SNPs genotyped in >66% of indivs; mean read depth for genotype call = 5..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.vcf --max-missing 0.66 --maf 0.05 --min-meanDP  5 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP05maxDP1000maf05
echo -e "\nKeep SNPs genotyped in >66% of indivs; mean read depth for genotype call = 10..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.vcf --max-missing 0.66 --maf 0.05 --min-meanDP 10 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP10maxDP1000maf05
echo -e "\nKeep SNPs genotyped in >66% of indivs; mean read depth for genotype call = 20..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.vcf --max-missing 0.66 --maf 0.05 --min-meanDP 20 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP20maxDP1000maf05
#the same as three above but with filtering (step F4) - finally used this... (Walter - usually used min-meanDP 10)
echo -e "\nKeep SNPs genotyped in >66% of indivs (from F4); mean read depth for genotype call = 5..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf --max-missing 0.66 --maf 0.05 --min-meanDP  5 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05
echo -e "\nKeep SNPs genotyped in >66% of indivs (from F4); mean read depth for genotype call = 10..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf --max-missing 0.66 --maf 0.05 --min-meanDP 10 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05
echo -e "\nKeep SNPs genotyped in >66% of indivs (from F4); mean read depth for genotype call = 20..."
vcftools --gzvcf "$species".SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf --max-missing 0.66 --maf 0.05 --min-meanDP 20 --max-meanDP 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05
#keep only one SNP per contig; however, no choice on which one, likely the first
echo -e "\nKeep only one SNP per contig (min-meanDP 5)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.meanDP05maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP05maxDP1000maf05thin
echo -e "\nKeep only one SNP per contig (min-meanDP 10)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.meanDP10maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP10maxDP1000maf05thin
echo -e "\nKeep only one SNP per contig (min-meanDP 20)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.meanDP20maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.meanDP20maxDP1000maf05thin
echo -e "\nKeep only one SNP per contig in F4 (min-meanDP 5)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05thin
echo -e "\nKeep only one SNP per contig in F4 (min-meanDP 10)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin
echo -e "\nKeep only one SNP per contig in F4 (min-meanDP 20)..."
vcftools --gzvcf "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05.recode.vcf --thin 1000 --recode --recode-INFO-all --out "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05thin

#copy data back
echo -e "\nCopying data back..."
cp "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05.recode.vcf ${DATADIR}/filtered.vcf
cp "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05thin.recode.vcf ${DATADIR}/filtered.vcf
cp "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf ${DATADIR}/filtered.vcf
cp "$species".SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05thin.recode.vcf ${DATADIR}/filtered.vcf
cp TotalRawSNPs.vcf ${DATADIR}/"$species".TotalRawSNPs.vcf
cp "$species".SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf ${DATADIR}/filtered.vcf

cp lowDP75.indv ${DATADIR}/filtered.vcf/"$species".lowDP75.indv
cp out.imiss ${DATADIR}/filtered.vcf/"$species".out.imiss
cp het_excess15.loci ${DATADIR}/filtered.vcf/"$species".het_excess15.loci
cp "$species"*.vcf ${DATADIR}/filtered.vcf
#cp dDocent_main.LOG ${DATADIR}/filtered.vcf/"$species".dDocent_main.LOG
#cp namelist ${DATADIR}/filtered.vcf/"$species".namelist
#cp reference.fasta ${DATADIR}/filtered.vcf/"$species".reference.fasta

#extract number of SNPs and number of contigs for all VCF files
echo -e "\nExtracting number of SNPs and number of contigs for all VCF files..."
files=(TotalRawSNPs.vcf
    ${species}.TotalRawSNPsPrim.vcf
    ${species}.SNPs.g5mac3.recode.vcf
    ${species}.SNPs.g5mac3minDP3bi.recode.vcf
    ${species}.SNPs.g5mac3minDP3birl15.recode.vcf
    ${species}.SNPs.g5mac3minDP3birl15rm75.recode.vcf
    ${species}.SNPs.g5mac3minDP3birl15rm75.recode.F1.vcf
    ${species}.SNPs.g5mac3minDP3birl15rm75.recode.F2.vcf
    ${species}.SNPs.g5mac3minDP3birl15rm75.recode.F3.vcf
    ${species}.SNPs.g5mac3minDP3birl15rm75.recode.F4.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05.recode.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05.recode.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05.recode.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP05maxDP1000maf05thin.recode.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP10maxDP1000maf05thin.recode.vcf
    ${species}.SNPs.g66mac3minDP3birl15rm75.F4.meanDP20maxDP1000maf05thin.recode.vcf
    )

for file in ${files[*]}; do
	echo "Logging loci and stacks in ${file}."
	echo "${file} $(grep -Eo '^dDocent\w*' ${file} | wc -l) $(grep -Eo '^dDocent\w*' ${file} | uniq | wc -l)" >> vcf.files_N.loci.stacks.txt
done

echo -e "\nCopying loci/contig statistics back..."
cp vcf.files_N.loci.stacks.txt ${DATADIR}/filtered.vcf/${species}.vcf.files_N.loci.stacks.txt

echo -e "\nFinished SNP filtering & statistics..."


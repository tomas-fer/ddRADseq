# ddRADseq data processing on [MetaCentrum](https://metavo.metacentrum.cz/en/index.html) using [dDocent](https://ddocent.com/)  
## 1. demultiplexing + trimming (script [`demultiplex_stacksMETA.sh`](demultiplex_stacksMETA.sh))  
using command `process_radtags` from [STACKS](https://catchenlab.life.illinois.edu/stacks/) (example of barcode [file](barcodesExample.txt))  
`process_radtags -1 ${file1} -2 ${file2} -b ${barcodes} -o demultiplexed --renz_1 ecoRI --renz_2 mspI -c –q`  
and [cutadapt](https://github.com/marcelm/cutadapt)  
`cutadapt --cores=$TORQUE_RESC_TOTAL_PROCS --cut 3 -U 5 --output "cutadapt/${fastq1}" --paired-output "cutadapt/${fastq2}" $fastq1 $fastq2`  

## 2. subsample to, e.g., 1,000,000 read pairs (script [`dDocentMETA_rawreadSubsampling.sh`](dDocentMETA_rawreadSubsampling.sh))  
using [seqtk](https://github.com/lh3/seqtk) command  
`seqtk sample -s100 ${file}.fq.gz ${subsam} > ${subsam}/${file}.fq`  

## 3. run dDocent (script [`dDocentMETA_fullRun.sh`](dDocentMETA_fullRun.sh))  
using command  
`dDocent config.file`  
see [config.file](config.file)  

## 4. SNP filtering (script [`SNPfilter_dDocentMETA.sh`](SNPfilter_dDocentMETA.sh))  
on `TotalRawSNPs.vcf` using [`vcfallelicprimitives`](https://github.com/vcflib/vcflib), [`vcftools`](https://vcftools.github.io/index.html) and [`vcffilter`](https://github.com/vcflib/vcflib)  

## 5. running sNMF ([LEA](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12382)) as implemented in [SambaR](https://github.com/mennodejong1986/SambaR)  
create PED/MAP and RAW/BIM files from filtered VCF using [`vcftools`](https://vcftools.github.io/index.html) and [`plink`](https://www.cog-genomics.org/plink/)  
`name="{your own prefix}"`  
`vcftools --gzvcf ${name}.vcf.gz --plink --out ${name}`  
`#correct PED file - 1st column should contain pop name (taken from the full name 'popname_samplename' in the second column)`  
`cut -f1 ${name}.ped | cut -f1 -d'_' > pops`  
`cut -f2- ${name}.ped > rest`  
`paste pops rest > ${name}.ped`  
`rm pops rest`  
`#calculate depth`  
`vcftools --gzvcf ${name}.vcf.gz --depth`  
`vcftools --gzvcf ${name}.vcf.gz --site-mean-depth`  
`#Convert PED/MAP to RAW/BIM`  
`plink --file ${name} --chr-set 95 --allow-extra-chr --make-bed --recode A --out ${name}`  

Run in R:  
`source("SAMBAR_v1.10.txt")`  
`getpackages()`  
`name="{your own prefix}"`  
`importdata(inputprefix=name)`  
`#filter data, i.e. here retain all SNPs (already filtered)`  
`filterdata(indmiss=1, snpmiss=1, min_mac=0, dohefilter=FALSE, snpdepthfilter=FALSE, min_spacing=0)`  
`findstructure(onlyLEA=TRUE, Kmax=7)`  

## 6. draw results as geographic maps with piecharts  
(scripts [`LEAmakePieMaps.sh`](LEAmakePieMaps.sh) and [`plotLEA_maps.R`](plotLEA_maps.R))  






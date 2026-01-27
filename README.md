# ddRADseq data processing on MetaCentrum using [dDocent](https://ddocent.com/)  
## 1. demultiplexing + trimming (script [`demultiplex_stacksMETA.sh`](demultiplex_stacksMETA.sh))  
using command `process_radtags` from [STACKS](https://catchenlab.life.illinois.edu/stacks/)  
`process_radtags -1 ${file1} -2 ${file2} -b ${barcodes} -o demultiplexed --renz_1 ecoRI --renz_2 mspI -c –q`  
and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
`java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 6 -phred33 demultiplexed/${file}.fq.gz trimmed/${file}.trimmed.fq.gz HEADCROP:10`  

## 2. subsample to, e.g., 1,000,000 read pairs (script [`dDocentMETA_rawreadSubsampling.sh`](dDocentMETA_rawreadSubsampling.sh))  
using [seqtk](https://github.com/lh3/seqtk) command  
`seqtk sample -s100 ${file}.fq.gz ${subsam} > ${subsam}/${file}.fq`  

## 3. run dDocent (script [`dDocentMETA_fullRun.sh`](dDocentMETA_fullRun.sh))  
using command  
`dDocent config.file`  
see [config.file](config.file)  

## 4. SNP filtering (script [`SNPfilter_dDocentMETA.sh`](SNPfilter_dDocentMETA.sh))  
on `TotalRawSNPs.vcf` using [`vcfallelicprimitives`](https://github.com/vcflib/vcflib), [`vcftools`](https://vcftools.github.io/index.html) and [`vcffilter`](https://github.com/vcflib/vcflib)  

## 5. running sNMF (LEA) as implemented in [SambaR](https://github.com/mennodejong1986/SambaR)  
create PED/MAP and RAW/BIM files using [`vcftools`](https://vcftools.github.io/index.html) and [`plink`](https://www.cog-genomics.org/plink/)  
Run in R:  
`source("SAMBAR_v1.10.txt")`  
`getpackages()`  
`importdata(inputprefix=name)`  
`#filter data, i.e. here retain all SNPs (already filtered)`  
`filterdata(indmiss=1,snpmiss=1,min_mac=0,dohefilter=FALSE,snpdepthfilter=FALSE,min_spacing=0)`  
`findstructure(onlyLEA=TRUE,Kmax=7)`  

## 6. draw results as geographic maps with piecharts  
(scripts [`LEAmakePieMaps.sh`](LEAmakePieMaps.sh) and [`plotLEA_maps.R`](plotLEA_maps.R))  






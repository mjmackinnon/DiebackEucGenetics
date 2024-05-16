#!/bin/bash

#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l jobfs=5GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

## run it like:
# qsub -v keeplist=lowmiss_geo_allspp Code/vcf_filter2.sh

module use /g/data/if89/apps/modulefiles
module load vcftools

echo ${keeplist}.keeplist

cd /scratch/xe2/jb5097/tmp/

vcftools --bcf Chr10.bcf \
	--keep ${keeplist}.keeplist \
	--remove-indels \
	--min-alleles 2 --max-alleles 2 \
	--max-missing 0.8 \
	--maf 0.01 --thin 200 \
	--min-meanDP 5 \
	--max-meanDP 25 \
	--recode --stdout | gzip -c > ${keeplist}_Chr10_snps2a_mis80_maf01_mDP5-18_thin200.vcf.gz

file=${keeplist}_Chr10_snps2a_mis80_maf01_mDP5-18_thin200

vcftools --gzvcf $file.vcf.gz --depth --out $file
vcftools --gzvcf $file.vcf.gz --site-mean-depth --out $file
vcftools --gzvcf $file.vcf.gz --het --out $file
vcftools --gzvcf $file.vcf.gz --relatedness --out $file
vcftools --gzvcf $file.vcf.gz --relatedness2 --out $file
vcftools --gzvcf $file.vcf.gz --missing-indv --out $file
#vcftools --gzvcf $file.vcf.gz --hardy --out $file
#vcftools --gzvcf $file.vcf.gz --TajimaD --out $file

echo "finished ${file}"



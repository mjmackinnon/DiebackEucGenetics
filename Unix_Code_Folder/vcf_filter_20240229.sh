#!/bin/bash
#PBS -l ncpus=24
#PBS -l mem=164GB
#PBS -l jobfs=24GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=04:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

## run it like:
# qsub Code/vcf_filter_20240229.sh

module use /g/data/if89/apps/modulefiles
module load vcftools

cd /scratch/xe2/jb5097/tmp/

in=/g/data/xe2/projects/paneuc2/final-variants/2024-02-29/mpileup~bwa~Emelliodora_sf2~melsider~filtered-default.vcf.gz

file=20240229_Chr10

vcftools --gzvcf $in \
	--chr Chr10 \
	--remove-indels \
	--min-alleles 2 --max-alleles 2 \
	--max-missing 0.8 \
	--min-meanDP 5 \
        --max-meanDP 25 \
	--maf 0.01 --thin 200 \
	--recode --stdout | gzip -c > ${file}.vcf.gz
# docs: https://vcftools.sourceforge.net/man_latest.html

vcftools --gzvcf $file.vcf.gz --depth --out $file
vcftools --gzvcf $file.vcf.gz --site-mean-depth --out $file
vcftools --gzvcf $file.vcf.gz --het --out $file
#vcftools --gzvcf $file.vcf.gz --relatedness --out $file
#vcftools --gzvcf $file.vcf.gz --relatedness2 --out $file
vcftools --gzvcf $file.vcf.gz --missing-indv --out $file
#vcftools --gzvcf $file.vcf.gz --hardy --out $file
#vcftools --gzvcf $file.vcf.gz --TajimaD --out $file

echo "finished ${file}"



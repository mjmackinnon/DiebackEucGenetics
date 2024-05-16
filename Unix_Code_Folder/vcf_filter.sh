#!/bin/bash

#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l jobfs=20GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

module use /g/data/if89/apps/modulefiles
module load vcftools

cd /scratch/xe2/jb5097/tmp/

vcftools --bcf Chr10.bcf \
	--remove-indels \
	--max-missing 0.8 \
	--maf 0.01 --thin 200 \
	--min-meanDP 5 \
	--max-meanDP 18 \
	--recode --stdout | gzip -c > Chr10_snps_mis80_maf01_mDP5-18_thin200.vcf.gz
	


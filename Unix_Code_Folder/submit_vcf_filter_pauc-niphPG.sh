for i in 01 02 03 04 05 06 07 08 09 10 11
do
echo Chr${i}
qsub -v keeplist=pauc_niph_PG,CHROM=Chr${i} Code/vcf_filter_pauc-niphPG.sh
done

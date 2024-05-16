## Fst scans all chromosomes

library(plyr)
library(ggplot2)

scans = ldply(list.files(pattern = "_biasnps_mis80_maf01_mDP5-25_fst_niph_pauc.weir.fst"), read.delim, header=TRUE)

table(scans$CHROM)

p = ggplot(scans, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 0.1) + 
  facet_wrap(~CHROM, nrow = 1, scales = 'free_x') + 
  labs(x = "Chr", y = "FST")
p


p <- ggplot(data = scans, aes(x = CHROM, y = WEIR_AND_COCKERHAM_FST)) + 
  geom_point(data = scans, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, fill = CHROM)) + 
  labs(x = "Chr", y = "FST")

p



ggplot(data = scans, aes(x = CHROM, y = WEIR_AND_COCKERHAM_FST)) + 
  geom_point(data = scans, 
             size = 0.1,
             aes(x = POS, y = WEIR_AND_COCKERHAM_FST, fill = CHROM)) + 
  labs(x = "Chr", y = "FST")


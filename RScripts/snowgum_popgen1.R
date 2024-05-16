## Lookin at population genetic structure in snowgums and related species using a sample of bilallelic SNPs from Chr 10

.libPaths( c( "/g/data/xe2/John/R/x86_64-pc-linux-gnu-library/4.3" , .libPaths() ) )

setwd("/scratch/xe2/jb5097/tmp/")

library(dplyr)
library(ggplot2)
library(stringr)

######
# PART 1 -- check data, get species names

file="filt2_A_Chr10_snps_mis80_maf01_mDP5-18_thin200"
# also try next:
# file = "filt2_A_Chr10_snps2a__mis80_maf01_mDP5-18_thin200"

path=paste0("/scratch/xe2/jb5097/tmp/")

# metadata:
metadata = read.csv2("snowgum_metadata.csv",
                     sep = ",", skip = 1)

## Check Fis distribution
het = read.delim2(paste0(path,file,".het"))
Fis = as.numeric(het$F)
hist(Fis, breaks = 20)

## missingness
imiss = read.delim2(paste0(path,file,".imiss"))
fmiss = as.numeric(imiss$F_MISS)
hist(fmiss, breaks = 20)

## sequence depth averaged across individuals
idepth = read.delim2(paste0(path,file,".idepth"))
idepth = as.numeric(idepth$MEAN_DEPTH)
hist(idepth, main = 'sequence depth averaged across individuals')  

## sequence depth averaged across sites
ldepth = read.delim2(paste0(path,file,".ldepth.mean")) %>% 
  mutate_at(c('POS', 'MEAN_DEPTH', 'VAR_DEPTH'), as.numeric)

summary(ldepth)
plot(ldepth$POS, ldepth$MEAN_DEPTH)
hist(ldepth$MEAN_DEPTH, main = 'sequence depth averaged across sites')

### join the sequencing summary data, merge it with sample metadata by sample name

## test whether names are unique
names = het$INDV
MONnames = names[grep("MON", names)]
length(MONnames) == length(unique(MONnames)) # not all MON names are unique!! technical reps?

BRYnames = names[grep("BRY", names)]
length(BRYnames) == length(unique(BRYnames)) # all BRY names are unique, phew. 

## A few steps to make sample names compatible between genetic data and metadata:

names_uniq = het$INDV
# for all the MON samples, only keep the first 8 characters. 
# This makes it compatible with metadata
sample = het$INDV
sample[grep("MON", sample)] = str_extract(sample[grep("MON", sample)], "^.{7}")
sample[grep("BRY", sample)] = str_sub(sample[grep("BRY", sample)], 1,-6)

sample_group = str_extract(names_uniq, "^.{3}")

df = as.data.frame(x = cbind(names_uniq, sample, sample_group, Fis, fmiss, idepth)) %>%
  mutate_at(c('Fis', 'fmiss', 'idepth'), as.numeric)

head(df)

plot(df$fmiss, df$idepth, col = factor(df$sample_group))
# Legend
legend("topright",
       legend = levels(factor(df$sample_group)),
       pch = 19,
       col = factor(levels(factor(df$sample_group))))

### Join the metadata with the genetic data for the samples that have metadata
# note: not all gen samps are in metadata
# note: some gen samps appear to be technical reps, these should get identical metadata

df2 <- df %>%
  left_join(metadata, by = 'sample')

# consolidate the species names
df2$Species2 <- if_else(df2$Species == "Eucalyptus pauciflora subsp niphophila", "niph",
                        if_else(df2$Species == "Eucalyptus pauciflora subsp pauciflora", "pauc",
                                if_else(df2$Species == "E. pauciflora", "pauc",
                                        if_else(df2$Species == "E. niphophila", "niph",
                                                if_else(df2$Species == "Eucalptus pauciflora hybrid pauc x niph", "paucXniph",
                                                        if_else(df2$Species == "E. radiata", "radi",
                                                                if_else(df2$Species == "E. stellulata", "stel",
                                                                        if_else(df2$Species == "E. delegatensis", "dele", "PROBLEM"))))))))

unique(df2$Species2)

write.table(df2, paste0(file, "_df2.tsv"), quote = F, row.names = F, sep = "\t")

plot(df2$fmiss, df2$idepth, col = factor(df2$Species2))
# Legend
legend("topright",
       legend = levels(factor(df2$Species2)),
       pch = 19,
       col = factor(levels(factor(df2$Species2))))

p = ggplot(df2, aes(x = Species2, y = Fis)) +
  geom_boxplot() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

## Plot a map of the samples and color by species

library(sf)
library(stars)

# read in the DEM:
DEMfile = paste0(path,"dem-9s.tif")
x = read_stars(DEMfile) # https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/66006
plot(x, axes = TRUE)

gdf = df2 %>%
  st_as_sf(coords=c("Longitude","Latitude"), crs = 4326) %>%
  st_transform(st_crs(x))

st_crs(gdf) == st_crs(x)

buf = 0.1
xmin = min(st_coordinates(gdf)[,1]) - buf
xmax = max(st_coordinates(gdf)[,1]) + buf
ymin = min(st_coordinates(gdf)[,2]) - buf
ymax = max(st_coordinates(gdf)[,2]) + buf
ROI = st_bbox(c(xmin = xmin, xmax = xmax, ymax = ymax, ymin = ymin), crs = st_crs(gdf))
ROI
xROI = st_crop(x, ROI, crop = TRUE)

plot(xROI, axes = TRUE)

ggplot() +
  geom_stars(data = xROI,
             downsample = 3) +
  theme_void() +
  scale_fill_viridis_c() +
  coord_equal() +
  geom_sf(data = gdf, aes(color = Species2))

# Add a column for geographic cluster
gdf$Lat <- st_coordinates(gdf)[,2]
gdf$Lon <- st_coordinates(gdf)[,1]
hist(gdf$Lat, breaks = 50)
b1 = -35
b2 = -32.5
abline(v = c(b1, b2), lty = 2)

gdf$geocluster <- if_else(gdf$Lat < b1, "South", 
                     if_else(gdf$Lat > b2, "North", "Middle"))

table(gdf$geocluster, gdf$Species2)


#######
## PART TWO -- Check genetic structure between recognised 'species'
# PCA, Fst, Admixture proportions

library(vcfR)

snps <- read.vcfR(paste0(file,".vcf.gz"),
                 checkFile = F,
                 #nrows = 5000,
                 convertNA = T) # don't check file because "bcf" in file header throws it off... 

# extract genotype data:
snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = F,
                             as.numeric = T,
                             convertNA = T,
                             return.alleles = F)

# rotate:
snps_num_t <- t(snps_num) 
snps_num_df <- data.frame(snps_num_t) 

# # check that sample names from VCF and from gdf are same:
# l1 = c("A", "B", "C")
# l2 = c("A", "Bx", "C")
# l1 == l2

row.names(snps_num_t) == gdf$names_uniq

length(gdf$Species2) == length(row.names(snps_num_t))

# check the genotype matrix
dim(snps_num_t)

# remove any columns(loci) with zero variance:
drop_loci = which(apply(snps_num_t, 2, var)==0)
length(drop_loci)

if(length(drop_loci)>0){
  message("Dropping some invariant loci now. Compare dims of snp matrix before and after...")
  snps_num_t2 <- snps_num_t[,-drop_loci]
  dim(snps_num_t2)
  SNPs_scaled <- scale(snps_num_t2)
} else {
  message("No invariant loci found. Now scaling")
  SNPs_scaled <- scale(snps_num_t)
}

pca_scaled <- prcomp(SNPs_scaled)

screeplot(pca_scaled, 
          ylab  = "Relative importance",
          main = "SNP Scree Plot")

prop.var <- summary(pca_scaled)$importance[2,]

# check:
pca_scaled$x[1:10,1:4]

pca_scores <- vegan::scores(pca_scaled)
pca_scores2 <- data.frame(gdf$Species2, 
                          gdf$geocluster,
                          gdf$names_uniq,
                          pca_scores)
colnames(pca_scores2)[1:2] <- c("Species", "Geocluster")


#pca_scores2[1:10,1:10]

dimx <- 1
dimy <- 2
ggpubr::ggscatter(data = pca_scores2,
                  y = paste0("PC", dimy),
                  x = paste0("PC", dimx),
                  color = "Species",
                  shape = "Geocluster",
                  alpha = 0.5,
                  xlab = paste0("PC", dimx, " (", round(prop.var[dimx]*100, 2), "% variation)"),
                  ylab = paste0("PC", dimy, " (", round(prop.var[dimy]*100, 2), "% variation)"))

dimx <- 3
dimy <- 4
ggpubr::ggscatter(data = pca_scores2,
                  y = paste0("PC", dimy),
                  x = paste0("PC", dimx),
                  color = "Species",
                  shape = "Geocluster",
                  alpha = 0.5,
                  xlab = paste0("PC", dimx, " (", round(prop.var[dimx]*100, 2), "% variation)"),
                  ylab = paste0("PC", dimy, " (", round(prop.var[dimy]*100, 2), "% variation)"))

# FST

# # Get sample lists, then pass it through VCFtools
# spp = "dele"
# geo = "South"
# maxN = 20
# 
# outlist <- gdf %>%
#   filter(geocluster == geo) %>%
#   filter(Species2 == spp) %>%
#   arrange(desc(idepth)) %>%
#   distinct(sample, .keep_all = FALSE) %>%
#   head(n = 20)
# 
# 
# subsample_pops <- function(gdf, spp, geo, maxN){
#   outlist <- gdf %>%
#     filter(geocluster == geo) %>%
#     filter(Species2 == spp) %>%
#     arrange(desc(idepth)) %>%
#     distinct(sample, .keep_all = TRUE) %>%
#     head(n = maxN)
#   
#   N = length(outlist$names_uniq)
#   
#   outfile = paste0(path, spp, "_", geo, "_n", N, ".popfile")
#   write(outlist$names_uniq, outfile)
# }
# 
# for(i in c("niph", "pauc","radi","stel","dele")){
#   subsample_pops(gdf, i, "South", 50)
# }
# 
# for(j in c("South","Middle", "North")){
#   subsample_pops(gdf, "pauc", j, 50)
# }

# open all the Fst files, get mean Fst, plot scan. 
fstfiles = list.files(pattern = "weir.fst")
# fix this

fstfile = fstfiles[10]
popA = strsplit(fstfile, '[.]')[[1]][1]
popB = strsplit(fstfile, '[.]')[[1]][2]
fst = read.delim2(fstfile)
FST = as.numeric(fst$WEIR_AND_COCKERHAM_FST)
Pos = as.numeric(fst$POS)
fst_mean = round(mean(FST, na.rm = TRUE), 3)
fst_median = round(median(FST, na.rm = TRUE), 3)

cbind(Pos, FST) %>%
  ggplot(aes(Pos/1000000, FST)) + 
  geom_point(alpha = 0.3) +
  ggtitle(paste0(popA, " X ", popB, ". mean:",fst_mean, " median: ", fst_median)) +
  xlab("Position (Mbp)")

# plot in windowed average
t <- as.data.frame(cbind(Pos, FST))
t$Pos_Mbp <- round(Pos/1000,0)

t %>%
  group_by(Pos_Mbp) %>%
  summarise(Fst = mean(FST)) %>%
  ggplot(aes(Pos_Mbp, Fst)) + 
  geom_point(alpha = 0.3) +
  ggtitle(paste0(popA, " X ", popB, ". mean:",fst_mean, " median: ", fst_median))

hist(FST, breaks = 100)

########
# PART THREE -- spatial genetic structure within 'species'

library(Rmisc)

# i   get equal matrices of pairwise genetic and geographic distance. 
# ii  heirachical clustering
# iii plot mean geog vs gen distance in distance bins (use code from Handroanthus analysis)

rel <- read.table(paste0(file,".relatedness2"), header = T)

get_geodist <- function(gdf){
  # generate geographic distance matrix for all genetic samples in a gdf
  d_gen <- as.matrix(st_distance(st_as_sf(gdf), st_as_sf(gdf)))
  rownames(d_gen) <- gdf$names_uniq
  colnames(d_gen) <- gdf$names_uniq
  d_gen[upper.tri(d_gen)] <- NA
  diag(d_gen) <- NA
  return(d_gen)
}

get_gendist <- function(rel, samps){
  
  rel_ <- rel %>%
    filter(INDV1 %in% samps) %>%
    filter(INDV2 %in% samps)
  
  rel_mat <- matrix(data = rel_$RELATEDNESS_PHI, 
                    nrow = length(unique(rel_$INDV1)), 
                    ncol = length(unique(rel_$INDV1)), 
                    byrow = TRUE,
                    dimnames = NULL)
  rownames(rel_mat) <- unique(rel_$INDV1)
  colnames(rel_mat) <- unique(rel_$INDV1)
  
  rel_mat[upper.tri(rel_mat)] <- NA
  diag(rel_mat) <- NA
  
  return(rel_mat)
}

unique(gdf$Species2)
unique(gdf$geocluster)

table(gdf$Species2, gdf$geocluster)

spp = "stel"
geoclust = "South"

gdf_ <- gdf %>%
  filter(Species2 == spp) %>%
  filter(geocluster == geoclust)

geo_mat <- get_geodist(gdf_)
geo_mat[1:10,1:10]
hist(geo_mat, breaks = 100)

gen_mat <- get_gendist(rel, gdf_$names_uniq)
gen_mat[1:10,1:10]
hist(gen_mat, breaks = 100, xlim = c(0,0.5))

colnames(geo_mat) == colnames(gen_mat)

# combine distances to df:
# Next, join the matrices into a datarame, remove the zeros, arrange by distance:
my.results <- data.frame(as.vector(geo_mat), as.vector(gen_mat))
names(my.results) <- c("Distance", "Phi")
my.results <- my.results %>%
  arrange(Distance) %>%
  na.omit()
head(my.results)

# take a quick look at relatedness decay:
plot(my.results$Distance, my.results$Phi)

# Now, bin the data:

bin_size= 100

meanPhi = mean(my.results$Phi)

dat_binned <- my.results %>%
  dplyr::group_by(dist_bin = cut(Distance, breaks = seq(0, max(Distance), bin_size))) %>%
  filter(!is.na(dist_bin)) %>%
  dplyr::summarise(Phi_mean = mean(Phi),
                   Phi_95CIupper = CI(Phi, ci = 0.95)[1],
                   Phi_95CIlower = CI(Phi, ci = 0.95)[3],
                   Phi_SD = sd(Phi),
                   N_pairs = n())
dat_binned$MidPoint <- seq(bin_size/2, (bin_size*nrow(dat_binned)) - bin_size, length.out = nrow(dat_binned))

dat_binned$MidPoint <- seq(from = bin_size/2, 
                           by = bin_size, length.out = nrow(dat_binned))

marks = c(0, 1,10,100,1000,10000,100000)
plot(dat_binned$MidPoint, dat_binned$Phi_mean,
     log = 'x',
     pch=20, col=rgb(0,0,0,alpha=0.9),
     xlab = "Pairwise Distance between individuals (m)",
     ylab = "Pairwise Relatedness (Phi)",
     main = paste0("Species: ",spp, "| group: ", geoclust, " | bin size: ", bin_size),
     ylim = c(0, 0.5),
     xaxt = "n"
     )
axis(1, at = marks, labels=format(marks, scientific=TRUE))

lines(dat_binned$MidPoint, dat_binned$Phi_mean)
lines(dat_binned$MidPoint, dat_binned$Phi_95CIlower, lty = 2)
lines(dat_binned$MidPoint, dat_binned$Phi_95CIupper, lty = 2)
abline(h = meanPhi, col = "red")

# ## plot with raw:
# 
plot(my.results$Distance, my.results$Phi,
     pch=20, col=rgb(0,0,0,alpha=0.05),
     xlab = "Pairwise Distance between individuals (m)",
     ylab = "Pairwise Relatedness (Phi)",
     main = paste0("Species: ",spp, " | bin size: ", bin_size))
#points(dat_binned$MidPoint, dat_binned$Phi_mean, pch=20, col=rgb(1,0,0,alpha=1))
lines(dat_binned$MidPoint, dat_binned$Phi_mean, lwd = 2, col = "red")
lines(dat_binned$MidPoint, dat_binned$Phi_95CIlower, lty = 2, col = "red")
lines(dat_binned$MidPoint, dat_binned$Phi_95CIupper, lty = 2, col = "red")
abline(h = meanPhi)
abline(h = 0, lty = 2)

#######
### Dendrogram

get_gendist2 <- function(rel, samps){
  
  rel_ <- rel %>%
    filter(INDV1 %in% samps) %>%
    filter(INDV2 %in% samps)
  
  rel_mat <- matrix(data = rel_$RELATEDNESS_PHI, 
                    nrow = length(unique(rel_$INDV1)), 
                    ncol = length(unique(rel_$INDV1)), 
                    byrow = TRUE,
                    dimnames = NULL)
  rownames(rel_mat) <- unique(rel_$INDV1)
  colnames(rel_mat) <- unique(rel_$INDV1)
  
  return(rel_mat)
}

library(vegan)
gen_dist_mat <- (get_gendist2(rel, gdf$names_uniq) - 0.5) * -1
gen_dist_mat <- as.dist(gen_dist_mat)
gen_dist_mat[1:10,1:10]
clust <- hclust(gen_dist_mat, method = 'average')

plot(clust, las = 1, 
     xlab="Sample", 
     ylab="Euclidean distance")

fit <- cascadeKM(gen_dist_mat, 1, 10, iter = 5000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)

### DENDROGRAM 2:
######
######

## Get some nice colors:
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)



library(ape)
dd <- dist(SNPs_scaled, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
hcp <- as.phylo(hc)

# plot dendrogram with nodes color coded by species, and labelled by species
hcp2 <- hcp
#hcp2$tip.label <- species.order
hcp2$tip.label <- paste(pca_scores2$Species, pca_scores2$Geocluster, pca_scores2$gdf.names_uniq, sep = "_")

# create a vector with one color per species name in the order of tips
colors <- paste(pca_scores2$Species, pca_scores2$Geocluster, sep = "_")
v2 <- c25[1:length(unique(colors))]

for(i in seq(1:length(unique(colors)))){
  sp <- unique(colors)[i]
  color <- v2[i]
  print(sp)
  print(color)
  colors[which(colors == sp)] <- color
}

pdf(file = paste0(file,"dendro.pdf"),
    width = 20, height = 30)
plot(hcp2, 
     tip.color = colors,
     label.offset = 0,
     cex = 0.5)
dev.off()

pdf(file = paste0(file,"dendro_unrooted.pdf"),
    width = 20, height = 20)
plot(hcp2, 
     tip.color = colors,
     label.offset = 0,
     cex = 0.5,
     type = 'unrooted')
dev.off()


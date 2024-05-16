## Lookin at population genetic structure in snowgums and related species using a sample of bilallelic SNPs from Chr 10

# ALL SPECIES 

setwd("/scratch/xe2/jb5097/tmp/")

library(dplyr)
library(ggplot2)
library(stringr)

######
# PART 1 -- check data, get species names

file="lowmiss_geo_allspp_Chr10_snps2a_mis80_maf01_mDP5-18_thin200"

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
  left_join(metadata, by = 'sample') %>%
  filter(Latitude != "")

# consolidate the species names
# df2$Species2 <- if_else(df2$Species == "Eucalyptus pauciflora subsp niphophila", "niph",
#                         if_else(df2$Species == "Eucalyptus pauciflora subsp pauciflora", "pauc",
#                                 if_else(df2$Species == "E. pauciflora", "pauc",
#                                         if_else(df2$Species == "E. niphophila", "niph",
#                                                 if_else(df2$Species == "Eucalptus pauciflora hybrid pauc x niph", "paucXniph",
#                                                         if_else(df2$Species == "E. radiata", "radi",
#                                                                 if_else(df2$Species == "E. stellulata", "stel",
#                                                                         if_else(df2$Species == "E. delegatensis", "dele", 
#                 if_else(df2$Species == "E. elata", "elat",
#                         if_else(df2$Species == "E. viminalis", "vimi",
#                                 if_else(df2$Species == "E. dives", "dive",
#                                         if_else(df2$Species == "E. rubida", "rubi",
#                                                 if_else(df2$Species == "E. fastigata", "fast",
#                                                         if_else(df2$Species == "E. dalrympleana", "dalr",
#                                                                 if_else(df2$Species == "E. aggregata", "aggr",
#                                                                         if_else(df2$Species == "E. canobolensis", "cano",
#                         "PROBLEM")))))))) ## not finished

df2$Species2 <- if_else(df2$Species == "Eucalyptus pauciflora subsp niphophila", "E. p. niphophila",
                        if_else(df2$Species == "Eucalyptus pauciflora subsp pauciflora", "E. pauciflora",
                                if_else(df2$Species == "Eucalptus pauciflora hybrid pauc x niph", "pauc X niph",
                                        if_else(df2$Species == "E. niphophila", "E. p. niphophila",
                                        df2$Species))))

table(df2$Species2)

df2 %>%
  filter(! Species2 %in% c("E. aggregata", "E. elata", "E. fastigata")) %>%
  ggplot(aes(x=fmiss, y = idepth, colour = Species2)) +
  geom_point(alpha = 0.5)

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

table(gdf$geocluster, gdf$Species)

#######
## PART TWO -- Check genetic structure between recognised 'species'
# PCA, Fst, Admixture proportions

library(vcfR)

snps <- read.vcfR(paste0(file,".vcf.gz"),
                  checkFile = F,
                  #nrows = 5000,
                  convertNA = T) # don't check file because "bcf" in file header throws it off... 

sample_nos <- match(gdf$names_uniq, df2$names_uniq)
snps <- snps[,c(1, sample_nos+1)] # add 1 so it doesnt take the format column!!


# extract genotype data:
snps_num <- vcfR::extract.gt(snps, 
                             element = "GT",
                             IDtoRowNames  = F,
                             as.numeric = T,
                             convertNA = F,
                             return.alleles = F)

# rotate:
snps_num_t <- t(snps_num) 
snps_num_df <- data.frame(snps_num_t) 

# check that sample names from VCF and from gdf are same:
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
                          pca_scores)
colnames(pca_scores2)[1:2] <- c("Species", "Geocluster")


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



# check species counts:
as.data.frame(table(pca_scores2$Species)) %>%
  arrange(desc(Freq))

# choose a species list for plotting
sp_list <- pca_scores2$Species

sp_list <- c("E. p. niphophila", "E. dives", "E. radiata", "E. aggregata", "E. stellulata",
             "E. blakelyi", "E. bridgesiana", "E. rubida",
             "E. viminalis", "E. delegatensis", "E. dalrympleana")

sp_list <- c("E. p. niphophila",  "E. stellulata",
             "E. blakelyi", "E. bridgesiana", "E. rubida",
             "E. viminalis", "E. dalrympleana")

pca_scores3 <- pca_scores2 %>%
  filter(Species %in% sp_list)

dimx <- 1
dimy <- 2
ggpubr::ggscatter(data = pca_scores3,
                  y = paste0("PC", dimy),
                  x = paste0("PC", dimx),
                  color = "Species",
                  #shape = "Geocluster",
                  alpha = 0.8,
                  palette = c25[1:length(sp_list)],
                  xlab = paste0("PC", dimx, " (", round(prop.var[dimx]*100, 2), "% variation)"),
                  ylab = paste0("PC", dimy, " (", round(prop.var[dimy]*100, 2), "% variation)"))

dimx <- 3
dimy <- 4
ggpubr::ggscatter(data = pca_scores3,
                  y = paste0("PC", dimy),
                  x = paste0("PC", dimx),
                  color = "Species",
                  #shape = "Geocluster",
                  alpha = 0.8,
                  palette = c25[1:length(sp_list)],
                  xlab = paste0("PC", dimx, " (", round(prop.var[dimx]*100, 2), "% variation)"),
                  ylab = paste0("PC", dimy, " (", round(prop.var[dimy]*100, 2), "% variation)"))

dimx <- 5
dimy <- 6
ggpubr::ggscatter(data = pca_scores3,
                  y = paste0("PC", dimy),
                  x = paste0("PC", dimx),
                  color = "Species",
                  #shape = "Geocluster",
                  alpha = 0.8,
                  palette = c25[1:length(sp_list)],
                  xlab = paste0("PC", dimx, " (", round(prop.var[dimx]*100, 2), "% variation)"),
                  ylab = paste0("PC", dimy, " (", round(prop.var[dimy]*100, 2), "% variation)"))

#### Important QUestion. How much of this pattern explained by shared missingness among species?

#### DENDROGRAM

library(ape)

# Compute distances and hierarchical clustering
dd <- dist(SNPs_scaled, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

plot(as.phylo(hc), cex = 0.6, label.offset = 0.5)

hcd <- as.dendrogram(hc)
plot(hcd, type = "rectangle", ylab = "Height")

# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# Customized plot; remove labels
plot(hcd, ylab = "Height", nodePar = nodePar, leaflab = "none")
# Horizontal plot
plot(hcd,  xlab = "Height",
     nodePar = nodePar, horiz = TRUE)

# the vector of species names in the order of dendrogram tips
species.order <- pca_scores2$Species[hc$order]

# create a vector with one color per species name in the order of tips
colors <- species.order
v2 <- c25[1:length(unique(pca_scores2$Species))]

for(i in seq(1:length(unique(species.order)))){
  sp <- unique(species.order)[i]
  color <- v2[i]
  print(sp)
  print(colour)
  colors[which(colors == sp)] <- color
}

# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.1, col = colors)
# Customized plot; remove labels
plot(hcd, ylab = "Height", nodePar = nodePar, leaflab = "none")
# Horizontal plot
plot(hcd, horiz = TRUE)

plot(dendrapply(hcd, colors))



library(ape)
plot(as.phylo(hc), cex = 0.9, label.offset = 1)
plot(as.phylo(hc), type = "unrooted")
clus5 = cutree(hc, 5)

mypal = c25[1:length(unique(pca_scores2$Species))]

hcp <- as.phylo(hc)
#hcp$tip.label <- species.order
plot(hcp, tip.color = mypal[clus5], label.offset = 1)

hcp



######
######
library(ape)
dd <- dist(SNPs_scaled, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
hcp <- as.phylo(hc)

# base plot with ape:
plot(hcp, label.offset = 1)

# get a vector of species names in the same order as dendrogram tips
species.order <- pca_scores2$Species[hc$order]

# # create a vector with one color per species name in the order of tips
# colors <- species.order
# v2 <- c25[1:length(unique(pca_scores2$Species))]
# 
# for(i in seq(1:length(unique(species.order)))){
#   sp <- unique(species.order)[i]
#   color <- v2[i]
#   print(sp)
#   print(color)
#   colors[which(colors == sp)] <- color
# }

# plot dendrogram with nodes color coded by species, and labelled by species
hcp2 <- hcp
#hcp2$tip.label <- species.order
hcp2$tip.label <- pca_scores2$Species

# create a vector with one color per species name in the order of tips
colors <- hcp2$tip.label
v2 <- c25[1:length(unique(colors))]

for(i in seq(1:length(unique(colors)))){
  sp <- unique(species.order)[i]
  color <- v2[i]
  print(sp)
  print(color)
  colors[which(colors == sp)] <- color
}

pdf(file = paste0(file,"dendro.pdf"),
    width = 20, height = 20)
plot(hcp2, 
     tip.color = colors,
     label.offset = 0,
     cex = 0.1)
dev.off()

pdf(file = paste0(file,"dendro_unrooted.pdf"),
    width = 20, height = 20)
plot(hcp2, 
     tip.color = colors,
     label.offset = 0,
     cex = 0.1,
     type = 'unrooted')
dev.off()


## Found that plotting the dendrogram with ape has it's own ordering convention, presuambly using hc$order.
## Therefore, give tip labels exactly corresponding to the snps rows, and colour vector with that. 
## possibly see nodelabels() to explore alternatives. 






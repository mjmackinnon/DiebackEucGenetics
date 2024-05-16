## JTB
# March 27 2024
# Merge physiology trait data with genetic sample meta data for snowgums collectde by Cal
# Approx. 30 E. p. niphophila and 30 E. p. pauciflora

P_meta <- read.csv2("/home/106/jb5097/Projects/SnowgumGenetics/Data/snowgum_phys_trait_data.csv",
                    sep = ",", skip = 1)
G_meta <- read.csv2("/home/106/jb5097/Projects/SnowgumGenetics/Data/snowgum_metadata.csv",
                    sep = ",", skip = 1) %>%
  dplyr::filter(grepl("BRY",sample))



new = str_sub(G_meta$sample, start = 5L, end = 7L)
new2 = str_replace_all(new, '[^[:alnum:]]', "")
G_meta$site.rep <- new


merged <- P_meta %>%
  left_join(G_meta, by = 'site.rep') %>%
  select(site.rep,Subspecies,leaf.size,sample,Species) %>%
  filter(!is.na(sample))

merged$leaf.size <- as.numeric(merged$leaf.size)

ggplot(merged, aes(x=Subspecies, y=leaf.size)) + 
  geom_boxplot()


## These are the samples to keep in the VCF. Next, match these up and filter based on sequencing depth/etc. 
merged$sample

## Read in "df2" which has output of vcftools and includes the genetic sample names. 
## 

d <- read.delim("filt2_A_Chr10_snps_mis80_maf01_mDP5-18_thin200_df2.tsv", sep = "\t") %>%
  dplyr::filter(grepl("BRY",sample))
d$site.rep <- sapply(strsplit(d$sample, "_"), "[[", 2)

merged2 = merged %>%
  left_join(d, by = 'site.rep') %>%
  filter(!is.na(Fis)) %>%
  select(site.rep, Subspecies, names_uniq, sample.y, Latitude, Longitude, idepth, fmiss)



### Map the distribution of samples:

library(sf)
library(stars)

# read in the DEM:
DEMfile = paste0(path,"dem-9s.tif")
x = read_stars(DEMfile) # https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/66006
plot(x, axes = TRUE)

gdf = merged2 %>%
  st_as_sf(coords=c("Longitude","Latitude"), crs = 4326) %>%
  st_transform(st_crs(x))

st_crs(gdf) == st_crs(x)

buf = 0.2
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
  geom_sf(data = gdf, aes(color = Subspecies))

# Sampling was done along a transect with multiple indiv in a small location.
# Seems that individual coordinates have high resolution
unique(df2$Latitude) 

## Save the samplelist

length(unique(merged2$site.rep)) # make sure no samples had two sequence libraries..
length(unique(merged2$names_uniq))

table(gdf$Subspecies)

PG_samples <- gdf$names_uniq

write(PG_samples, paste0(path,"pauc-niph_PG.keeplist"))

for(i in c("niph", "pauc", "paucXniph")){
  samps <- d$names_uniq[which(d$Species2==i)]
  print(i)
  print(length(samps))
  write(samps, paste0(path,i,".popfile"))
}






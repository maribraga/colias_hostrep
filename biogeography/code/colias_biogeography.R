# Colias Biogeography - code in Quarto document

library(tidyverse)
library(sf)
sf::sf_use_s2(FALSE)

bioregions <- st_read("./biogeography/data/bioregions.geojson") 

palette <- c("NT" = "#FF8C01",
             "NA" = "#DA3541",
             "EP" = "#00A2FF",
             "PTP" = "#F8C700",
             "WP" = "#013459",
             "AF" = "#EF9FBA",
             "SA" = "grey50")

ggplot(bioregions) +
  geom_sf(aes(fill = biome, col = biome)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme_bw()


## Sarah's workflow ----

library(rnaturalearth)
library(printr)
library(regionfeatures)

#data <- load_shapefile("~/repos/sswiston-regionfeatures/regionfeatures_package/data/shapefiles/sam_regions_final.shp", "Region")
#class(data)
#areas <- areas(data, filepath="~/repos/sswiston-regionfeatures/regionfeatures_package/data", format="CSV")

# comment to Sarah: it only works if you are in the parent directory to /data
setwd("/Users/mari/repos/colias_hostrep/biogeography/data/regionfeatures")
data_7regions <- load_shapefile("./data/shapefiles/bioregions.shp", "biome")

# comment to Sarah: have to change the projection to equal area
# Sarah uses an equal area projection for South America
# I'm using Gall's projection, which works better for the whole world
library(rgdal)
data_7regions_equal <- data_7regions
data_7regions_equal$data <- spTransform(data_7regions$data, CRS("+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# as sf object
data_7regions_sf <- st_as_sf(data_7regions$data)

# library(tmap)
# qtm(data_7regions_sf)

ggplot(data_7regions_sf) +
  geom_sf(aes(fill = biome, col = biome)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values = palette) +
  theme_bw()


#### test without west palearctic ##
# 
# bio_test <- bioregions[c(4,7),]
# bio_test <- st_transform(bio_test, crs = "+proj=cea +lon_0=0 +x_0=0 +y_0=0 +lat_ts=45 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# st_write(bio_test, dsn = "americas.shp")
# 
# test <- load_shapefile("americas.shp", "biome")
# data_7regions <- test


# Generate areas

areas <- areas(data_7regions_equal, filepath="./data", format="CSV")
# Contents
summary(areas)
# Quantitative areas
prmatrix(areas$areas, rowlab=rep("",1))
# Classification (1 = above mean, 0 = below mean)
prmatrix(areas$classification, rowlab=rep("",1))


# Generate annual mean temperatures (BIO1)

bio1 <- bioclim(data_7regions_equal, c(1), filepath="./data", format="CSV")
# comment to Sarah: it's a list of 1 element that contains the list of interest
bio1 <- bio1$bio1
# Contents
summary(bio1)
# Mean temperatures
prmatrix(bio1$mean, rowlab=rep("",1))
# Classification (1 = above mean, 0 = below mean)
prmatrix(bio1$classification, rowlab=rep("",1))
# Difference in mean temperatures (absolute value)
prmatrix(bio1$mean_diff)
# Sameness (1 = regions share classification, 0 = regions do not share classification)
prmatrix(bio1$sameness)


# Generate distances

# find bug
distances <- distances(data_7regions_equal, filepath="./data", format="CSV")
# Contents
summary(distances)
# Adjacency (1 = adjacent, 0 = non-adjacent)
prmatrix(distances$adjacency)
# Mean distances (mean distance between any two points in each region)
prmatrix(distances$mean)


# Generate altitudes

# find bug
alts <- altitudes(data_7regions_equal, filepath="./data", format="CSV")
# Contents
summary(alts)
# Mean altitudes
prmatrix(alts$mean, rowlab=rep("",1))
# Classification (1 = above mean, 0 = below mean)
prmatrix(alts$classification, rowlab=rep("",1))
# Difference in mean altitude (absolute value)
prmatrix(alts$mean_diff)
# Sameness (1 = regions share classification, 0 = regions do not share classification)
prmatrix(alts$sameness)




regions <- data_7regions_equal

# distances()
reorder = "NONE" 
filepath="./data"
format = "CSV"


# altitudes()
reorder = "NONE" 
filepath = "." 
format = "NONE"
n_regions <- regions$n_regions
region_names <- regions$region_names
suppress_ar = 1
output = list()
name <- "altitudes"
elevation_raster <- elevatr::get_elev_raster(regions$data, z = 1)
elev <- raster::projectRaster(from = elevation_raster, crs = CRS("+init=epsg:4326"))
                                             neg_to_na = TRUE, override_size_check = TRUE
descriptions <- read.csv("data/variable_descriptions.csv", 
                         as.is = TRUE)
desc <- descriptions[which(descriptions$Feature == name),]
description <- desc$Description
units <- desc$Units
raster = list()
raster = list(data = elevation_raster, name = name, description = description, 
              units = units, type = "quantitative", tag = "raster")

output <- analyze_raster(regions = regions, reorder = reorder, 
                         raster = raster, filepath = filepath, format = format, 
                         quiet = TRUE)
alts <- output

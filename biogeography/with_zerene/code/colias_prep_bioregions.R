# Get rasters of bioeographical regions

library(sf)
library(terra)
library(tmap)
library(tmaptools)
library(tidyverse)
library(rgdal)
library(broom)
library(geojsonio)
 
# The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
# which was just loaded, will retire in October 2023.
# Please refer to R-spatial evolution reports for details, especially
# https://r-spatial.org/r/2023/05/15/evolution4.html.


# Nicolas regions for animated plot
nico_biomes <- geojson_read("~/repos/colias_hostrep/biogeography/data/Nico_biomes.geojson",  what = "sp")
geojson_write(nico_biomes, file = "~/repos/colias_hostrep/biogeography/data/Nico_biomes.geojson")
#st_read("~/repos/colias_hostrep/biogeography/data/Nico_biomes.geojson")

# 'fortify' the data to get a dataframe format required by ggplot2
nico_biomes_f <- tidy(nico_biomes)

# Plot it
ggplot() +
  geom_polygon(data = nico_biomes_f, aes( x = long, y = lat, group = group, fill = id)) +
  geom_polygon(data = pth_f, aes( x = long, y = lat, group = group), fill="black", alpha = 0.6) +
  theme_void() +
  coord_map()

# Rasmus version
tmap_options(check.and.fix = TRUE)
qtm(nico_biomes)

## Geojson files
# Pan-Tibetan Highlands
qtp <- st_read("/Users/mari/repos/colias_hostrep/biogeography/data/Pan-Tibetan-Highlands/GeoJson/pan_tibetan_highlands.geojson")
qtm(qtp)

# Pan-Tibetan Highlands + adjacent mountains
tp <- st_read("/Users/mari/repos/colias_hostrep/biogeography/data/Pan-Tibetan-Highlands/GeoJson/tibetan_plateau_adjacent.geojson")
qtm(tp)


## Shape files
# Pan-Tibetan Highlands

pth <- readOGR( 
  dsn = paste0(getwd(),"/biogeography/data/Pan-Tibetan/PTH/Shapefile/") , 
  layer = "Pan-Tibetan_Highlands_P",
  verbose = FALSE
)
# 'fortify' the data to get a dataframe format required by ggplot2
pth_f <- tidy(pth)

# Pan-Tibetan Highlands + adjacent mountains
ptha <- readOGR( 
  dsn = paste0(getwd(),"/biogeography/data/Pan-Tibetan/PTH_expanded/Shapefile/") , 
  layer = "PTHadj_P",
  verbose = FALSE
)
# 'fortify' the data to get a dataframe format required by ggplot2
ptha_f <- tidy(ptha)


## background for plot

world <- ne_countries(scale = "medium", returnclass = "sf")
eurasia <- subset(world, continent %in% c("Europe", "Asia"))

# not working
# box <- c(xmin=-30,xmax=170,ymin=-85,ymax=85)
# eurasia <- st_crop(eurasia, box)


# Plot it
ggplot() +
  geom_sf(data = eurasia) +
  #geom_polygon(data = ptha_f, aes( x = long, y = lat, group = group), fill="darkorange", color="darkorange") +
  geom_polygon(data = pth_f, aes( x = long, y = lat, group = group), fill="#69b3a2", color="#69b3a2", alpha = 0.6) +
  coord_sf(
    xlim = c(-30, 170),
    ylim = c(-10, 90)
  ) +
  theme_bw() 

### Combine Nicolas regions and PTP:
## Merge 8-9, 3-6, drop 1 and 5, add PTP
## Move border between 4 and 7 to west end of PTP (keeping both regions contiguous)

# Rasmus lekstuga

sf::sf_use_s2(FALSE)

nico_biomes <- geojson_read("./biogeography/data/my_biomes copy.geojson",  what = "sp") %>% 
  st_as_sf()

# works with geojson_read
ggplot() +
  geom_polygon(data = tidy(nico_biomes), aes( x = long, y = lat, group = group, fill = id)) +
  theme_void() +
  coord_map()

all <- st_read("./biogeography/data/my_biomes.geojson") 
Nam <- st_read("./biogeography/data/NA.geojson")
Sea <- st_read("./biogeography/data/SEA.geojson") 
ptp <- st_read("./biogeography/data/Pan-Tibetan/PTH/GeoJson/pan_tibetan_highlands.geojson") %>% 
  st_union() %>% 
  st_as_sf() %>% 
  rename(geometry = x) %>% 
  mutate(biome = "PTP")
division <- st_read("./biogeography/data/division.shp") 

EP <- all %>% filter(biome =="EP") 
west_east_asia <- st_intersection(EP, division)
qtm(newEP)

WP <- all %>% filter(biome =="WP") 

newWP <- st_union(WP, west_east_asia) %>% 
  st_difference(.,ptp)
qtm(newWP)

newEP <- st_difference(EP, west_east_asia)
newEP <- st_difference(newEP, ptp)
qtm(newEP)

bioregions <- bind_rows(Sea, newEP, newWP, Nam, ptp,
                        filter(all, biome %in% c("NT","AF")))
bioregions <- bioregions %>% 
  select(biome, geometry)

bioregions_s <- st_simplify(bioregions, preserveTopology = TRUE)
qtm(bioregions_s, fill = "biome")

#st_write(bioregions, "./biogeography/data/bioregions.geojson")
#st_write(bioregions, "./biogeography/data/bioregions.shp")

palette <- c()

ggplot(bioregions) +
  geom_sf(aes(fill = biome)) +
  scale_fill_manual(values = palette) +
  theme_bw()

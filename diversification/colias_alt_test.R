library(here)


## Tree

tree <- read.tree(here("data_prepare/data/processed/Colias_full_tree_clean.tre"))

## Elevation
elev_original <- read.csv(paste0("./data_prepare/data/original/Colias_elevation.csv"))

elev_wide <- elev_original %>% 
  separate(elevation, into = c("min", "max"), "-") %>% 
  mutate(name = str_remove(species,"_MSF\\d+"),
         min = as.numeric(min),
         max = case_when(is.na(max) ~ min, TRUE ~ as.numeric(max)))

# duplicates <-  elev_wide %>% 
#   filter(str_detect(name, "\\d$")) %>% 
#   pull(name)

elevation <- elev_wide %>%
  #mutate(name = str_remove(name, "\\d$")) %>% 
  select(name, min, max) %>% 
  distinct() %>% 
  mutate(range = max-min)

#elevation <- filter(elevation, name != "Colias_interior_")

elevation_tree <- elevation %>% 
  filter(name %in% tree$tip.label) %>% 
  mutate(name = factor(name, levels = tree$tip.label))
  
## Geo

original_states <- read.csv(here("biogeography/without_zerene/data/colias_geo_tbl.csv")) %>% 
  arrange(factor(species, levels = tree_geo$tip.label)) %>% 
  rename(name = species, distribution = range)

data_el_dist <- left_join(elevation_tree, original_states) %>% 
  mutate(name = factor(name, levels = tree$tip.label))


## plot

palette <- c("#EF9FBA", "#6E4672", "#00D164", "#AA74A1", "#DA3541", "#FF8C01", 
             "#009F52", "#00A2FF", "#347FD5", "#013459", "#F8C700", "#097C65")

p <- ggtree(tree) + theme_tree2() #+ geom_tiplab(size = 2) + xlim(0,3.8)

pcol <- p %<+% data_el_dist + 
  geom_tippoint(aes(color = distribution)) +
  scale_color_manual(values = palette) +
  theme(legend.position = "none")

elev_plot <- ggplot(data_el_dist) +
  geom_linerange(aes(y = name, 
                     xmin = min, 
                     xmax = max, 
                     col = distribution)) +
  scale_color_manual(values = palette) +
  theme_bw() +
  labs(y = "", x = "Elevation (m)")

revts(pcol) + elev_plot


## Es-sim test

min <- pull(elevation_tree, min)
names(min) <- elevation_tree$name

max <- pull(elevation_tree, max)
names(max) <- elevation_tree$name

range <- pull(elevation_tree, range)
names(range) <- elevation_tree$name

source("./diversification/essim.r")

essim_min <- essim(tree, min, return.es = TRUE)
essim_max <- essim(tree, max)
essim_range <- essim(tree, range)

essim_min$rho
essim_min$`P Value`
essim_max
essim_range

par(mfrow = c(3,1))
plot(min, essim_min$es)
plot(max, essim_min$es)
plot(range, essim_min$es)



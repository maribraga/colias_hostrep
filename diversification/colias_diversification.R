#'---
#'title: "Colias diversification"
#'author: "Mariana P Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#+ setup, include = FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

#'-------------
#'
#' Script for analyses performed in ... 
#' *...*
#'

#' ## Set up
#' For this script we'll need  

#+ packages
library(evolnets)
library(MCMCpack)
library(coda)
library(kdensity)
library(igraph)
library(bipartite)
library(ape)
library(treeio)
library(ggtree)
library(patchwork)
library(tidyverse)
library(purrr)
  
#' ## 
#'

#' ** **
#' 
# read files
path_data <- "./diversification/data/"

tree <- read_tree_from_revbayes(paste0(path_data,"tree_final_Rev.tre"))
plot(tree, show.node.label = T)

## Elevation
elev_original <- read.csv(paste0(path_data, "colias_elevation.csv"))

elev_wide <- elev_original %>% 
  separate(elevation, into = c("min", "max"), "-") %>% 
  mutate(name = str_remove(species,"_MSF\\d+"),
         min = as.numeric(min),
         max = as.numeric(max))

duplicates <-  elev_wide %>% 
  filter(str_detect(name, "\\d$")) %>% 
  pull(name)

elevation <- elev_wide %>%
  mutate(name = str_remove(name, "\\d$")) %>% 
  select(name, min, max) %>% 
  distinct() %>% 
  mutate(range = max-min)

ggplot(elevation) +
  geom_density(aes(min), col = "grey50") +
  geom_density(aes(max), col = "grey10") +
  geom_density(aes(range), col ="blue")


# plot elevational range
min_order <- elevation  %>% 
  arrange(min) %>% 
  pull(name)

elevation %>% 
  ggplot() +
  geom_linerange(aes(x = factor(name, levels = min_order), ymin = min, ymax = max))

elevation_tree <- elevation %>% 
  filter(name %in% tree$tip.label) %>% 
  mutate(name = factor(name, levels = tree$tip.label))
  
elevation_tree %>% 
  ggplot() +
  geom_linerange(aes(x = factor(name, levels = tree$tip.label), ymin = min, ymax = max))

ggplot(elevation_tree) +
  geom_density(aes(min), col = "grey50") +
  geom_density(aes(max), col = "grey10")

elev_dropped <- setdiff(elevation$name, elevation_tree$name)


## Geographic distribution

# subjective division
distribution <- read.csv(paste0(path_data, "colias_distribution.csv")) %>% 
  mutate(name = str_remove(species,"_MSF\\d+"),
         name = str_remove(name, "\\d$")) %>%
  select(name, areas) %>% 
  distinct()
  
distribution_tree <- distribution %>% 
  filter(name %in% tree$tip.label) %>%
  mutate(name = factor(name, levels = tree$tip.label))
  
dist_dropped <- setdiff(distribution$name, distribution_tree$name)

intersect(elev_dropped, dist_dropped)



## Combine data

data_el_dist <- left_join(elevation_tree, distribution_tree)

p <- ggtree(tree) + theme_tree2() #+ geom_tiplab() + xlim(0,3.5)
pcol <- p %<+% data_el_dist + 
  geom_tippoint(aes(color=areas)) +
  theme(legend.position = "none")

elev_plot <- ggplot(data_el_dist) +
  geom_linerange(aes(y = name, 
                     xmin = min, 
                     xmax = max, 
                     col = areas)) +
  theme_bw() +
  labs(y = "", x = "Elevation (m)")

revts(pcol) + elev_plot


# how many species per area?
data_el_dist %>% 
  group_by(areas) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n))


## Summarize distribution by realms

neot <- c("A", "B")
near <- c("C", "D")
eastp <- c("E", "G", "H")
westp <- c("F")
sea <- c("I")
africa <- c("J")
# combinations: eastp and westp, eastp and westp and near, eastp and near, eastp and sea
east_asia <- c(eastp, sea)
palearctic <- c(eastp, westp)
nam_asia <- c(near, eastp)
arctic <- c(near, eastp, westp)

all_states <- c("neot","near","eastp","westp","sea","africa","east_asia","palearctic","nam_asia","arctic")

for(n in all_states){
  full <- cross(rep(list(get(n)), length(get(n)))) %>% 
    simplify_all() %>% 
    lapply(sort) %>% 
    lapply(unique) %>% 
    lapply(paste, collapse = '') %>% 
    unique() %>% 
    unlist()
  
  assign(paste0(n,"_full"), full)
}


data_el_dist_nico <- data_el_dist %>% 
  mutate(distribution = case_when(
    # single regions
    areas %in% neot_full ~ "Neotropical",
    areas %in% near_full ~ "Nearctic",
    areas %in% eastp_full ~ "Palearctic_east",
    areas %in% westp_full ~ "Palearctic_west",
    areas %in% sea_full ~ "Southeast_Asian",
    areas %in% africa_full ~ "African",
    # combinations
    areas %in% east_asia_full ~ "East_Asian",
    areas %in% palearctic_full ~ "Palearctic",
    areas %in% nam_asia_full ~ "N_American-Asian",
    areas %in% arctic_full ~ "Arctic",
    TRUE ~ NA_character_
    )
  )

# order ranges so that colors are more intuitive
# there is no species only in SE Asia
data_el_dist_nico_plot <- data_el_dist_nico %>% 
  mutate(Range = factor(distribution, levels = c("Neotropical","Nearctic","N_American-Asian","Arctic","Palearctic_east","East_Asian","Palearctic","Palearctic_west","African")))

palette <- c("#487D65","#DA3541","#AA74A1","#6E4672","#347FD5","#FCA63E","#007ea7","#003459","#91B2A3")

p_nico <- p %<+% data_el_dist_nico_plot + 
  geom_tippoint(aes(color = Range)) +
  scale_color_manual(values = palette) +
  theme(legend.position = "none")

elev_nico_plot <- ggplot(data_el_dist_nico_plot) +
  geom_linerange(aes(y = name, 
                     xmin = min, 
                     xmax = max, 
                     col = Range),
                 linewidth = 2) +
  scale_color_manual(values = palette) +
  theme_bw() +
  labs(y = "", x = "Elevation (m)")

p_nico + elev_nico_plot





# Define order of area combinations

library(ape)
library(tidyverse)
library(here)

tree <- read.tree(here("biogeography/without_zerene/data/tree_colias_geo.tre"))

range_tbl <- read.csv(here("biogeography/without_zerene/data/colias_geo_tbl.csv")) %>% 
  arrange(factor(species, levels = tree$tip.label))

range_matrix <- read.csv(here("biogeography/without_zerene/data/colias_geo_mtx.csv"), row.names = 1)

area_codes <- rep(1, ncol(range_matrix))
names(area_codes) <- colnames(range_matrix)
area_codes

state_labels <- read.csv(here("biogeography/without_zerene/server_ignore/output/run1/state_labels.txt"), colClasses = c("numeric","character")) %>% 
  separate_wider_position(range, widths = area_codes)

states <- state_labels %>% 
  pivot_longer(2:8, names_to = "area", values_to = "present") %>% 
  filter(present == 1) %>% 
  select(-present) %>% 
  group_by(state) %>% 
  reframe(combination = str_c(area, collapse = "+"))

sort(unique(states$combination))

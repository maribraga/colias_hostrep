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
library(phytools)
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
path_data <- "./data/"

treez <- read.tree(paste0(path_data,"full_tree_clean.tre"))
tree <- drop.tip(treez, "Zerene_cesonia")
plot(tree, show.node.label = T, cex = 0.5)
axisPhylo()

is.binary(tree)
is.ultrametric(tree) # why is it not ultrametric?
is.rooted(tree)


#### hide ----
tree_ultra <- force.ultrametric(tree)
plot(tree_ultra, show.node.label = T, cex = 0.5)
axisPhylo()

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
eastp <- c("E") 
qtp <- c("G", "H") # new
westp <- c("F")
sea <- c("I")
africa <- c("J")
# combinations: eastp and westp, eastp and westp and near, eastp and near, eastp and sea, west and qtp, east and qtp, east and west and qtp
east_asia <- c(eastp, sea, qtp)
palearctic_low <- c(eastp, westp)
palearctic <- c(eastp, westp, qtp)
east_qtp <- c(eastp, qtp) 
west_qtp <- c(westp, qtp)
nam_asia <- c(near, eastp)
arctic <- c(near, eastp, westp)


#all_states <- c("neot","near","eastp","westp","sea","africa","east_asia","palearctic","nam_asia","arctic")
all_states <- c("neot","near","eastp","qtp","sea","westp","africa","east_asia","palearctic_low","palearctic","east_qtp","west_qtp","nam_asia","arctic")

for(s in all_states){
  full <- cross(rep(list(get(s)), length(get(s)))) %>% 
    simplify_all() %>% 
    lapply(sort) %>% 
    lapply(unique) %>% 
    lapply(paste, collapse = '') %>% 
    unique() %>% 
    unlist()
  
  assign(paste0(s,"_full"), full)
}

# go from less inclusive to more inclusive
data_el_dist_nico <- data_el_dist %>% 
  mutate(distribution = case_when(
    # single regions
    areas %in% neot_full ~ "Neotropical",
    areas %in% near_full ~ "Nearctic",
    areas %in% eastp_full ~ "Palearctic_east",
    areas %in% qtp_full ~ "QTP",
    areas %in% sea_full ~ "Southeast_Asian",
    areas %in% westp_full ~ "Palearctic_west",
    areas %in% africa_full ~ "African",
    # combinations
    areas %in% east_qtp_full ~ "East_QTP",
    areas %in% west_qtp_full ~ "West_QTP",          
    areas %in% east_asia_full ~ "East_Asian",
    areas %in% palearctic_low_full ~ "Palearctic_notQTP",
    areas %in% palearctic_full ~ "Palearctic",
    areas %in% nam_asia_full ~ "N_American-Asian",
    areas %in% arctic_full ~ "Arctic",
    TRUE ~ NA_character_
    )
  )

## Make a matrix with areas and distribution color
# order ranges so that colors are more intuitive
# there is no species only in SE Asia or in East_QTP
spp_areas <- distribution %>% 
  mutate(areas = gsub("([A-Z])", "\\1 ", areas)) %>% 
  separate_wider_delim(areas, delim = " ", names_sep = "", too_few = "align_start") %>% 
  pivot_longer(!name, names_repair = "minimal")

data_el_dist_nico_plot <- data_el_dist_nico %>% 
  mutate(Range = factor(distribution, levels = c("Neotropical",
                                                 "Nearctic",
                                                 "N_American-Asian",
                                                 "Arctic",
                                                 "Palearctic_east",
                                                 "East_Asian",
                                                 "QTP",
                                                 "West_QTP",
                                                 "Palearctic",
                                                 "Palearctic_notQTP",
                                                 "Palearctic_west",
                                                 "African")))

palette <- c("#FF8C01",
             "#DA3541",
             "#AA74A1",
             "#6E4672",
             "#00A2FF",
             "#00D164",
             "#F8C700",
             "#097C65",
             "#009F52",
             "#347FD5",
             "#013459",
             "#EF9FBA")

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




#### RPanda ----

library(RPANDA)

#treeColias <- readRDS()

# time of MRCA
tot_time <- max(node.depth.edgelength(tree))

# sampling fraction - total of 66 species in the tree, but only 50 with host data
tot_described_species <- 79
sam_frac <- Ntip(tree)/tot_described_species 

#### _Rate shifts ----

missing <- c("Colias_harfordii", "Colias_occidentalis", "Colias_aegidii", "Colias_adelaidae", "Colias_felderi", "Colias_erschoffi", "Colias_sagartia", "Colias_aias", "Colias_leechi", "Colias_nilagiriensis", "Colias_marnoana", "Colias_ponteni", "Colias_dubia")
n_miss <- length(missing) 
#clades <- c(rep(1,7), rep(2,20), rep(3,17), rep(4,12), rep(5,11), rep(NA, n_miss))
ns <- c(rep("south",7), rep("north",60), rep("north", n_miss))
asia <- c(rep("americas",27),rep("asia",40), rep("americas",2),rep("asia",11))

groups <- data.frame(Species = c(tree$tip.label, missing), 
#                     clade = clades,
                     north_south = ns,
                     east_west = asia)

#### __quick start ----

library(ape)
library(RPANDA)

library(phytools)
library(ape)
library(phangorn)
library(parallel)
library(ParallelLogger)
library(RColorBrewer)
library(devtools)
library(pspline)

tree <- readRDS("tree_Colias.rds") 
groups <- readRDS("groups.rds")

sfrac <- get.sampling.fractions(tree, groups, plot = TRUE, cex = 0.3)
node_comb <- get.comb.shift(tree, groups, sfrac)

par(mfrow=c(2,2))
for(i in 1:4){
  plot.phylo.comb(phylo = tree, 
                data = groups, 
                sampling.fractions = sfrac, 
                combi = node_comb[i], 
                edge.width = 2, 
                label.offset = 0.3, 
                lad = F, 
                cex = 0.4)
  axisPhylo()
}


#### __crown shifts ----

shifts <- shift.estimates(phylo = tree,
                          data = groups,
                          sampling.fractions = sfrac,
                          comb.shift = node_comb,
                          models = c("BCST","BCST_DCST","BVAR",
                                    "BVAR_DCST","BCST_DVAR", "BVAR_DVAR"),
                          backbone.option = "crown.shift",
                          multi.backbone = "all",
                          Ncores = 4)

str(shifts, max.level = 1)
head(shifts$total)
shifts$backbones$`75/`

par(mfrow=c(1,1))
plot.phylo.comb(phylo = tree, 
                data = groups, 
                sampling.fractions = sfrac,
                shift.res = shifts,
                combi = 1, 
                edge.width = 2, 
                label.offset = 0.3, 
                lad = F, 
                cex = 0.4)
axisPhylo()

rates <- div.rates(phylo = tree,
                   shift.res = shifts,
                   combi = 1, part = "all")
rates
time <- -c(tot_time, floor(tot_time):0)

rates_by_clade <- tibble(clade = NA, time = NA, rate = NA)
for(i in 1:length(rates)){
  for(j in 1:length(time)){
    rates_by_clade[i,1] <- names(rates)[i]
    rates_by_clade[i,2] <- 
  }
}


#plot(time, rates[[length(rates)]][1,], type = "l", col = "blue", ylim = c(0, 1.3))
#lines(time, rates[[length(rates)]][2,], type = "l", col = "red")

plot(time, rates[[length(rates)]][1,], type = "l",
     #ylim = c(0,1.3),
     col = "blue",
     las = 1, xlab = "Time (Myrs)",
     ylab = "Rates (Events/ Lineages / Myr)",
     main = "Diversification rates of the backbone")
lines(time, rates[[length(rates)]][2,], type = "l", col = "red")
legend("topleft", legend = c("Speciation rate", "Extinction rate"),
       col = c("blue", "red"), lty = 1, bty = "n")

plot(time, rates[[length(rates)]][1,] - rates[[length(rates)]][2,], type = "l",
     ylab = "Diversification rate")
lines(time, rates[[1]][1,], type = "l", col = "red")
legend("topright", legend = c("Backbone", "North_clade"),
       col = c("black", "red"), lty = 1, bty = "n")

paleodiversity <- paleodiv(phylo = tree,
                                   data = groups,
                                   sampling.fractions = sfrac,
                                   shift.res = shifts)

plot(time, paleodiversity,
     bty = "n", las = 1, ylim = c(0, 100),
     main = "Paleodiversity dynamic of Colias",
     xlab = "Time (Myrs)", ylab = "Number of species", lwd = 2, type = "l", lty = 1)


#### __stem shifts ----

shifts_stem <- shift.estimates(phylo = tree,
                          data = groups,
                          sampling.fractions = sfrac,
                          comb.shift = node_comb,
                          models = c("BCST","BCST_DCST","BVAR",
                                     "BVAR_DCST","BCST_DVAR", "BVAR_DVAR"),
                          backbone.option = "stem.shift",
                          multi.backbone = "all",
                          Ncores = 4)

str(shifts_stem, max.level = 1)
head(shifts_stem$total)

plot.phylo.comb(phylo = tree, 
                data = groups, 
                sampling.fractions = sfrac,
                shift.res = shifts_stem,
                backbone.option = "stem.shift",
                combi = 1, 
                edge.width = 2, 
                label.offset = 0.3, 
                lad = F, 
                cex = 0.4)
axisPhylo()

rates_stem <- div.rates(phylo = tree,
                   shift.res = shifts_stem,
                   combi = 1, part = "all")
rates_stem
time <- -c(tot_time, floor(tot_time):0)

rates_by_clade <- tibble(clade = NA, time = NA, rate = NA)
for(i in 1:length(rates)){
  for(j in 1:length(time)){
    rates_by_clade[i,1] <- names(rates)[i]
    rates_by_clade[i,2] <- 
  }
}

plot(time, rates_stem[[length(rates_stem)]][1,], type = "l")
lines(time, rates_stem[[1]][1,], type = "l", col = "red")


#### __save objects to send to Helene ----

saveRDS(tree, "tree_Colias.rds")
saveRDS(groups, "groups.rds")
saveRDS(shifts, "shifts.rds")
saveRDS(shifts, "shifts_BD.rds")
saveRDS(shifts, "shifts_Bcst.rds")


#### __Compare results ----

shifts_B <- readRDS("shifts_B.rds")
shifts_BD <- readRDS("shifts_BD.rds")
shifts_Bcst <- readRDS("shifts_Bcst.rds")

rates_Bcst <- div.rates(phylo = tree,
                   shift.res = shifts_Bcst,
                   combi = 1, part = "all")
rates_Bvar <- div.rates(phylo = tree,
                   shift.res = shifts_Bvar,
                   combi = 1, part = "all")
rates_BD <- div.rates(phylo = tree,
                   shift.res = shifts_BD,
                   combi = 1, part = "all")
rates_Bcst
rates_BD


#### _Fit birth-death ----

# no extinction
f.lamb <- function(t,y){y[1]} # functional form for lambda
f.mu <- function(t,y){0} # functional form for mu, which is not being estimated here
lamb_par <- c(0.09)
mu_par <- c()
result_cst <- fit_bd(tree, tot_time, f.lamb, f.mu, lamb_par, mu_par,
                    f = sam_frac, cst.lamb=TRUE, fix.mu=TRUE)
result_cst$model <- "pure birth with constant speciation rate"
result_cst

# constant lambda and mu
f.lamb <- function(t,y){y[1]} # functional form for lambda
f.mu <- function(t,y){y[1]} # functional form for mu
lamb_par <- c(0.08)
mu_par <- c(0.001)
result_cst2 <- fit_bd(tree, tot_time, f.lamb, f.mu, lamb_par, mu_par,
                    f = sam_frac, cst.lamb=TRUE, cst.mu=TRUE)
result_cst2$model <- "constant speciation and extinction rates"
result_cst2

# it might output a negative value for the parameter, but just take the absolute value


# exponential lambda, no mu
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
result_exp <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     f=sam_frac,expo.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
result_exp$model <- "pure birth with exponential variation in speciation rate"
result_exp
plot_fit_bd(result_exp, tot_time)


# exponential lambda, constant mu
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.05, 0.01)
mu_par<-c(0.002)
result_exp_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     f=sam_frac,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
result_exp_cst$model <- "exponential variation in speciation rate and constant extinction"
result_exp_cst

# exponential mu, constant lambda
f.mu <-function(t,y){y[1] * exp(y[2] * t)}
f.lamb <-function(t,y){y[1]}
mu_par<-c(0.05, 0.01)
lamb_par<-c(0.002)
result_cst_exp <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     f=sam_frac,cst.lamb=TRUE,expo.mu=TRUE,dt=1e-3)
result_cst_exp$model <- "exponential variation in extinction rate and constant speciation"
result_cst_exp

AICc <- c(result_cst$aicc, result_cst2$aicc, result_exp$aicc, result_exp_cst$aicc, result_cst_exp$aicc)
names(AICc) <- c(result_cst$model, result_cst2$model, result_exp$model, result_exp_cst$model, result_cst_exp$model)

deltaAICc <- AICc - min(AICc)
deltaAICc

plot_fit_bd(result_cst_exp, tot_time)






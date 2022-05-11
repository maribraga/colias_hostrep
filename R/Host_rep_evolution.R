#'---
#'title: "Colias host repertoire evolution"
#'author: "Mariana P Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script for analyses performed in ... 
#' *...*
#'

#' ## Set up
#' For this script we'll need *evolnets* to analyze the output posterior distribution for character history.
#' You can install it from GitHub:

#+ eval = FALSE
# install.packages("devtools")
devtools::install_github("maribraga/evolnets")

#' We also need other packages

#+ message = FALSE

library(evolnets)
library(MCMCpack)
library(coda)
library(kdensity)
library(igraph)
library(ape)
library(treeio)
library(ggtree)
library(patchwork)
library(bipartite)
library(tidyverse)

#' # Parameter estimates
#' The first thing we need to do is to check that independent MCMC chains have converged to the same posterior distribution. Both chains have achieved an effective sample size ESS > 300 for every parameter.
#'

#' **log files**
#' 
# read files
log1 <- read.table("./R/data/out.4.2s.beta.colias.log", header = TRUE)
log2 <- read.table("./R/data/out.5.2s.beta.colias.log", header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")

#' **Convergence test**
#' 
gelman.diag(mcmc.list(as.mcmc(log1), as.mcmc(log2)))

#' **Parameter estimates**
#' Then, we can plot the posterior distributions. since both chains have samples from the same posterior distribution, we can use only one of the log files for that.
#' 
parameters <- log1 %>% 
  pivot_longer(cols = 2:5, names_to = "parameter", values_to = "value")

ggplot(parameters, aes(parameter, value)) +
  geom_violin() +
  stat_summary(fun.data = "median_hilow", color = "red", size = 0.5) +
  theme_bw()

parameters %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))

#' # Character history
#' Now we can move on to reconstruct the history of evolution of host repertoire across de Colias phylogeny.
#' 
#' ## Data
#' First we read in the phylogenetic trees for butterflies and plants.
#' Then, we read in the interaction matrix.
#'

#' **Trees**
#' We use `read_tree_from_revbayes` to read the Colias tree because this file was exported from RevBayes and contains the node labels given by RevBayes. This will be very important in the analysis!
tree <- read_tree_from_revbayes("./R/data/tree_Rev.tre")
host_tree <- read.tree("./R/data/host_tree.tre")

#' **Extant network**
#' This matrix contains 0s and 2s because in the host repertoire model in RevBayes, there are 3 possible states (0,1,2), where 1 means "potential host" and 2 means "actual host". We used the 2-state model for the reconstruction in RevBayes, so we are only interested in the 0s and 2s, no potential host.
#' 
matrix <- read.csv("./R/data/matrix.csv", row.names = 1) %>% as.matrix()





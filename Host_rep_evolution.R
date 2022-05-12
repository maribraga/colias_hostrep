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

#' ## Parameter estimates
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

#' All values are very close to 1, so we are good to go.
#' 

#' **Parameter estimates**
#' 
#' Now, we can plot the posterior distributions. Since both chains have samples from the same posterior distribution, we can use only one of the log files for that.
#' 
#+ densities, dpi = 300
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

#' ## Character history
#' Let's move on to reconstruction of the history of evolution of host repertoire across the Colias phylogeny.
#' 
#' ### Data
#'

#' **Trees**
#' 
#' We use `read_tree_from_revbayes` to read the Colias tree because this file was exported from RevBayes and contains the node labels given by RevBayes. This will be very important in the analysis!
tree <- read_tree_from_revbayes("./R/data/tree_Rev.tre")
host_tree <- read.tree("./R/data/host_tree.tre")

#' **Extant network**
#' 
#' This matrix contains 0s and 2s because in the host repertoire model in RevBayes, there are 3 possible states (0,1,2), where 1 means "potential host" and 2 means "actual host". We used the 2-state model for the reconstruction in RevBayes, so we are only interested in the 0s and 2s, no potential host.
#' 
matrix <- read.csv("./R/data/matrix.csv", row.names = 1) %>% as.matrix()


#' **Read in .history.txt files**
#' 
#' We'll use the *evolnets* function `read_history()` to read a file outputed from RevBayes with sampled histories during MCMC
#'  
history <- read_history("./R/data/out.4.2s.beta.colias.history.txt", burnin = 0.2)

#' ### Number of events and effective rate of evolution
#' 

# Estimated number of events across the Colias phylogeny
count_events(history)

# How many events were host gains and how many host losses?
count_gl(history)

# Considering the number of events and the total length of the Colias phylogeny, 
# what is the estimated rate of host repertoire evolution?
effective_rate(history, tree)


#' ### Extant and ancestral networks
#' 
#' **Time points of interest (ages)**
#' 
#' The first step to reconstruct ancestral networks is to define the time points during Colias diversification that we want to look at.
#' 

#+ tree, fig.width = 6, fig.height = 8, dpi = 300
# visually determine interesting time points to reconstruct ancestral networks
plot(tree, show.node.label = T, cex = 0.5)
axisPhylo()

# choose time points
ages <- c(0,0.4,0.8,1.2)

#' **Interaction probability at given ages**
#' 
#' Now we calculate the posterior probability of interaction between each host and each extant butterfly at each age in `ages`.
#' 
at_ages <- posterior_at_ages(history, ages, tree, host_tree)
pp_at_ages <- at_ages$post_states

#' **Summarize probabilities into networks**
#' 
# Make binary or weighted networks? Discard interactions with posterior probability < threshold.
# We chose to reconstruct weighted networks with a threshold of 

# ##### temporary fix ######
pp_at_ages <- lapply(pp_at_ages, abind::adrop, drop = 3)
# ###

summary_networks <- get_summary_network(pp_at_ages, pt = 0.8)

#+ eval = FALSE
# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree)
#  save R object to be able to use same module configuration in the future
saveRDS(modules_at_ages, "R/R_objects/modules_at_ages_pp80.rds")

#+ 
# using a saved object
modules_at_ages <- readRDS("R/R_objects/modules_at_ages_pp80.rds")


#' **Plot networks**
#' 

#+ ancestral_nets, fig.width = 15, fig.height = 13, warning = F, dpi = 300
# get the information needed from `modules_at_ages`
matched_modules <- modules_at_ages$matched_modules$nodes_and_modules_per_age

p_nets <- plot_ancestral_networks(summary_networks, matched_modules, tree) 
wrap_plots(p_nets$plot, nrow = 2)

#+ extant_graph, fig.width = 8, fig.height = 6.5, warning = F, dpi = 300
p_nets$plot[[4]]

#+ plotmoduleweb, fig.width = 7, fig.height = 7, warning = F, dpi = 300
mod_ext <- modules_at_ages$original_modules$moduleWeb_objects$`0`
plotModuleWeb(mod_ext, labsize = 0.4)






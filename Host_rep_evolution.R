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
#' The first thing we need to do is to check that independent MCMC chains have converged to the same posterior distribution. Both chains have achieved an effective sample size ESS > 200 for every parameter.
#'

#' **log files**
#' 
# read files
path_out<- "ignore/2s_new_tree/output/"
log1 <- read.table(paste0(path_out,"out.2.2b.colias.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.3.2b.colias.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")

#' **Convergence test**
#' 
its <- seq(50000,550000,1000)

chain1 <- filter(log1, iteration %in% its) %>% as.mcmc()
chain2 <- filter(log2, iteration %in% its) %>% as.mcmc()

gelman.diag(mcmc.list(chain1, chain2))

#' All values are very close to 1, so we are good to go.
#' 

#' **Parameter estimates**
#' 
#' Now, we can plot the posterior distributions. Since both chains have samples from the same posterior distribution, we can use only one of the log files for that.
#' 
#+ densities, out.width = '70%', dpi = 300
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

#' **Bayes factor**
#' 
#' The parameter called `beta` defines whether the phylogenetic distance between hosts affects the probability of gaining new hosts. When `beta = 0`, phylogenetic distances between hosts do not matter, hence all hosts are equally likely to be gained. We use Bayes factor to test whether the inferred `beta` is significantly different from 0. If the factor is < 1, it means that the phylogenetic distances do not matter during host gain events.
#'  
d_prior <- dexp(x=0, rate=1)

kd_beta <- kdensity(x = log1$beta, 
                    kernel='gamma', 
                    support=c(0,Inf), 
                    bw = 0.02)
max = kd_beta(0)

(BF <- d_prior/max)

/*
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

#+ tree, fig.width = 6, fig.height = 8, out.width = '80%', dpi = 300
# visually determine interesting time points to reconstruct ancestral networks
plot(tree, show.node.label = T, cex = 0.5)
axisPhylo()

# choose time points
ages <- c(0,0.4,0.8,1.2)

#' I've chosen 1.2, 0.8, and 0.4, which means 1.2 Ma, 800 and 400 thousand years ago. 
#' 
#' **Interaction probability at given ages**
#' 
#' Now we calculate the posterior probability of interaction between each host and each extant butterfly at each age in `ages`.
#' 
at_ages <- posterior_at_ages(history, ages, tree, host_tree)
pp_at_ages <- at_ages$post_states

#' **Summarize probabilities into ancestral networks**
#' 
#' Make binary or weighted networks? Discard interactions with posterior probability < threshold.
#' We chose to reconstruct weighted networks with a threshold of 
#' 

#+ echo = FALSE
# ##### temporary fix ######
pp_at_ages <- lapply(pp_at_ages, abind::adrop, drop = 3)
# ###

#+ 
summary_networks <- get_summary_network(pp_at_ages, pt = 0.8)

#+ eval = FALSE
# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree)
#  save R object to be able to use same module configuration in the future
saveRDS(modules_at_ages, "R/R_objects/modules_at_ages_pp80.rds")

#+ 
# using a saved object instead
modules_at_ages <- readRDS("R/R_objects/modules_at_ages_pp80.rds")


#' **Plot networks**
#' 

#+ ancestral_nets, fig.width = 15, fig.height = 13, dpi = 300, warning = F
# get the information needed from `modules_at_ages`
matched_modules <- modules_at_ages$matched_modules$nodes_and_modules_per_age

pal <- scales::hue_pal()(11)[c(11,2,5,8,9,1,10,3,4,6,7)]
  
p_nets <- plot_ancestral_networks(summary_networks, matched_modules, tree, palette = pal)
wrap_plots(p_nets$plot, nrow = 2)

/*
#+ extant_graph, fig.width = 8, fig.height = 6.5, dpi = 300, warning = F
p_nets$plot[[4]]

#+ plotmoduleweb, fig.width = 7, fig.height = 7, dpi = 300, warning = F
mod_ext <- modules_at_ages$original_modules$moduleWeb_objects$`0`
plotModuleWeb(mod_ext, labsize = 0.4)
*/


#' **Plot networks as matrices**
#' 
#' An alternative way to visualize the networks is by plotting them as matrices.
#' 

#+ ancestral_matrices, fig.width = 12, fig.height = 14, dpi = 300, warning = F
# get the modules for the extant network
mod_ext <- modules_at_ages$original_modules$moduleWeb_objects$`0`

# create plots with a matrix for each age
p_mod_matrix_ages <- list()

for(i in seq_along(ages)){
  
  a <- ages[i]
  
  if(a != 0){
    net <- summary_networks[[as.character(a)]] %>% as.matrix()
    
    mod_df <- modules_at_ages$matched_modules$nodes_and_modules_per_age %>% 
      filter(age == a) %>% 
      select(name, module, type) 
    
    plot <- plot_module_matrix(net, mod_df)
    
  } else{
    
    mod_list <- listModuleInformation(mod_ext)[[2]]
    host_mods <- lapply(mod_list, function(x) data.frame(host = x[[2]]))
    host_mods <- dplyr::bind_rows(host_mods, .id = 'host_module')
    mod_order <- host_mods$host
    
    para_mods <- lapply(mod_list, function(x) data.frame(parasite = x[[1]]))
    para_mods <- dplyr::bind_rows(para_mods, .id = 'parasite_module')
    mod_order_para <- para_mods$parasite
    
    plot <- plot_module_matrix(matrix, mod_ext, 
                       parasite_order = mod_order_para, 
                       host_order = mod_order)
  }
  
  p_mod_matrix_ages[[i]] <- plot
  
}

# define the layout of the plot
layout <- c(
  area(3,3,7,5),
  area(3,1,7,2),
  area(1,4,2,5),
  area(1,2)
)
/*
  plot(layout)
*/

# plot!
wrap_plots(p_mod_matrix_ages, ncol = 2, design = layout, guides = "keep")




/*  # hide this for now. not necessary and $plot is wrong
#' **Module validation**
#' 
#' We use modules to facilitate visualisation but sometimes they also have biological meaning. This step tells us whether the modules in the ancestral networks are robust across MCMC samples or an artifact of the summary networks.
#' 

samples_at_ages <- c(
  lapply(at_ages$samples[1:3], function(x) x[1:3,,]),
  at_ages$samples[4]
)
mod_samples <- modules_from_samples(samples_at_ages)

mod_val <- support_for_modules(mod_samples, modules_at_ages)
mod_val$plot
mod_val$pairwise_membership
mod_val$mean_support

*/

  
#' ### Ancestral states at internal nodes of Colias phylogeny
#' We also wanted to do a tradition ancestral state reconstruction (ASR), calculating  interaction probabilities at internal nodes of the Colias tree. And now that we have defined modules for the extant network, we can also use them to group hosts in the ASR. 
#' 
#' **Interaction probability at internal nodes**
#' 

#+ ancestral_states, fig.width = 15, fig.height = 13, dpi = 300, warning = F
at_nodes <- posterior_at_nodes(history, tree, host_tree)
p_asr <- plot_module_matrix2(matrix, at_nodes, tree, host_tree, modules = mod_ext, threshold = 0.9)
p_asr

*/

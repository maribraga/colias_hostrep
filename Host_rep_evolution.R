#'---
#'title: "Colias host repertoire evolution"
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
#' For this script we'll need *evolnets* to analyze the output posterior distribution for character history.
#' You can install it from GitHub:

#+ eval = FALSE
# install.packages("devtools")
devtools::install_github("maribraga/evolnets")

#' We also need other packages

#+ packages
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
log1 <- read.table(paste0(path_out,"out.5.2b.colias.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.2.2b.colias.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")

#' **Convergence test**
#' 
its <- seq(20000,200000,500)

chain1 <- filter(log1, iteration %in% its)
chain2 <- filter(log2, iteration %in% its)

gelman.diag(mcmc.list(as.mcmc(chain1), as.mcmc(chain2)))

effectiveSize(chain1)

#' All values are very close to 1, so we are good to go.
#' 

#' **Parameter estimates**
#' 
#' Now, we can plot the posterior distributions. Since both chains have samples from the same posterior distribution, we can use only one of the log files for that.
#' 
#+ densities, out.width = '50%', dpi = 300
parameters <- chain1 %>% 
  pivot_longer(cols = 2:5, names_to = "parameter", values_to = "value")

plot_param <- ggplot(parameters, aes(parameter, value)) +
  geom_violin(col = "grey40") +
  stat_summary(fun.data = "median_hilow", color = "#E76F51", size = 0.5) +
  theme_bw()

plot_param

parameters %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))

#' **Bayes factor**
#' 
#' The parameter called `beta` defines whether the phylogenetic distance between hosts affects the probability of gaining new hosts. When `beta = 0`, phylogenetic distances between hosts do not matter, hence all hosts are equally likely to be gained. We use Bayes factor to test whether the inferred `beta` is significantly different from 0. If the factor is < 1, it means that the phylogenetic distances do not matter during host gain events.
#'  
d_prior <- dexp(x=0, rate=1)

kd_beta <- kdensity(x = chain1$beta, 
                    kernel='gamma', 
                    support=c(0,Inf), 
                    bw = 0.02)
max = kd_beta(0)

(BF <- d_prior/max)


#' ## Character history
#' Let's move on to reconstruction of the history of evolution of host repertoire across the Colias phylogeny.
#' 
#' #### Data
#'

#' **Trees**
#' 
#' We use `read_tree_from_revbayes` to read the Colias tree because this file was exported from RevBayes and contains the node labels given by RevBayes. This will be very important in the analysis!
path_data <- "ignore/2s_new_tree/data/"
path_evol <- "ignore/2s_new_tree/evolnets/"

tree <- read_tree_from_revbayes(paste0(path_evol,"tree_final_Rev.tre"))
host_tree <- read.tree(paste0(path_data,"host_tree.tre"))

#' **Extant network**
#' 
#' This matrix contains 0s and 2s because in the host repertoire model in RevBayes, there are 3 possible states (0,1,2), where 1 means "potential host" and 2 means "actual host". We used the 2-state model for the reconstruction in RevBayes, so we are only interested in the 0s and 2s, no potential host.
#' 
matrix <- read.csv(paste0(path_data,"matrix_phylo_timetree.csv"), row.names = 1) %>% as.matrix()

#' **Read in .history.txt files**
#' 
#' We'll use the *evolnets* function `read_history()` to read a file outputed from RevBayes with sampled histories during MCMC
#'  
history <- read_history(paste0(path_out, "out.5.2b.colias.history.txt"), burnin = 0) %>% 
  filter(iteration %in% its)


#' #### Number of events and effective rate of evolution
#' 

# Estimated number of events across the Colias phylogeny
count_events(history)

# How many events were host gains and how many host losses?
gl <- count_gl(history)
# percentage of gains and losses
gl/sum(gl)

# Considering the number of events and the total length of the Colias phylogeny, 
# what is the estimated rate of host repertoire evolution?
effective_rate(history, tree)


#' #### Modules of the extant (present-day) network
#' 

#+ modules, eval = FALSE
# find modules
mod <- mycomputeModules(matrix)

#+ modules_rds, include = FALSE
mod <- readRDS(paste0(path_evol, "7mods.rds"))


#' ### Extant and ancestral networks
#' 
#' **Time points of interest (ages)**
#' 
#' The first step to reconstruct ancestral networks is to define the time points during Colias diversification that we want to look at.
#' 

#+ tree, fig.width = 6, fig.height = 8, out.width = '60%', dpi = 300
# visually determine interesting time points to reconstruct ancestral networks
pt <- ggtree(tree) + 
  geom_tiplab() + 
  geom_nodelab(size = 3, color = "grey40") + 
  theme_tree2()

pt_rev <- revts(pt) + xlim(c(-3, 0.8))

pt_rev_ages <- pt_rev + geom_vline(xintercept = c(-2.1,-1.4,-0.7,0), col = "blue")

# choose time points
ages <- c(0,0.7,1.4,2.1)

#' I've chosen 2.1, 1.4, and 0.7, which means 2.1 Ma, 1.4 Ma and 700 thousand years ago. 
#' 
#' **Interaction probability at given ages**
#' 
#' Now we calculate the posterior probability of interaction between each host and each extant butterfly at each age in `ages`.
#' 
at_ages <- posterior_at_ages(history, ages, tree, host_tree)

#' **Summarize probabilities into ancestral networks**
#' 
#' Make binary or weighted networks? Discard interactions with posterior probability < threshold.
#' We chose to reconstruct weighted networks with a threshold of 
#' 

summary_networks <- get_summary_networks(at_ages, threshold = 0.9)

#+ match_modules, eval = FALSE
# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree, extant_modules = mod)

#+ match_mod_read, echo  FALSE
#set.seed(19)
#saveRDS(modules_at_ages, paste0(path_evol,"matched_modules_90_9m_seed19.rds"))
modules_at_ages <- readRDS(paste0(path_evol,"matched_modules_90_9m_seed19.rds"))

# modules and submodules
mods_pal <- sort(unique(modules_at_ages$matched_modules$original_and_matched_module_names$module_name))
pal <- c("#2A9D8F","#8AB17D","#E9C46A",
          "#E89A5E","#E76F51","#D1495B",
          "#DF838F","#256C6A","#203A44")
pal_extant <- c("#2A9D8F","#8AB17D","#E9C46A",
                 "#D1495B","#DF838F","#256C6A","#203A44")


#' **Plot ancestral networks**
#' 

#+ ancestral_nets, fig.width = 12, fig.height = 10, dpi = 300
p_nets <- plot_ancestral_networks(summary_networks, modules_at_ages, tree, palette = pal, node_size = 2)
wrap_plots(p_nets, heights = c(2,3), guides = "collect")



#' **Plot ancestral networks as matrices**
#' 
#' An alternative way to visualize the networks is by plotting them as matrices.
#' 
  
#+ anc_nets_matrix, fig.width = 7, fig.height = 22, dpi = 300, out.width = '80%'
# create plots with a matrix for each age
p_mod_matrix_ages <- list()

for(i in seq_along(ages)){
  
  a <- ages[i]
  
  net <- summary_networks[[as.character(a)]] %>% as.matrix()
    
  mod_df <- modules_at_ages$matched_modules$nodes_and_modules_per_age %>% 
      filter(age == a) %>% 
      select(name, module, type)
    
  if(a != 0){
    
    plot <- plot_extant_matrix(net, mod_df) + 
      scale_fill_manual(values = pal, breaks = mods_pal) +
      labs(title = paste0(a," Ma"))
    
  } else{
    
    # names of matched modules for extant network
    match_mod_ext <- modules_at_ages$matched_modules$nodes_and_modules_per_age %>% 
      filter(age == 0) %>% 
      select(name, module, type)
    mod_order <- match_mod_ext %>% 
      filter(type == "host") %>% 
      arrange(module) %>% 
      pull(name)
    mod_order_para <- match_mod_ext %>% 
      filter(type == "symbiont") %>% 
      arrange(module) %>% 
      pull(name)
    
    p_ext_mods <- plot_extant_matrix(matrix, match_mod_ext, 
                                     parasite_order = mod_order_para, 
                                     host_order = mod_order)
    p_ext_mods + scale_fill_manual(values = pal_extant)
    
    plot <- plot_extant_matrix(net, mod_df, 
                               parasite_order = mod_order_para, 
                               host_order = mod_order) + 
      scale_fill_manual(values = pal, breaks = mods_pal) +
      labs(title = paste0(a," Ma"))
  }
  
  p_mod_matrix_ages[[i]] <- plot
  
}

# define the layout of the plot
layout <- c(
  area(1,1,4,4),
  area(5,1,7,3),
  area(8,1,9,2),
  area(10,1,10,2)
)

# plot!
wrap_plots(p_mod_matrix_ages, design = layout, guides = "keep")


#' **Module validation**
#' 
#' We use modules to facilitate visualisation but sometimes they also have biological meaning. This step tells us whether the modules in the ancestral networks are robust across MCMC samples or an artifact of the summary networks.
#' 
samples_at_ages <- get_sampled_networks(at_ages)

#+ mod_samp, eval = FALSE
mod_samples <- modules_from_samples(samples_at_ages)

#+ mod_samp_read, echo = FALSE
#saveRDS(mod_samples, paste0(path_evol, "mod_samples.rds"))
mod_samples <- readRDS(paste0(path_evol, "mod_samples.rds"))
# took 56 min!

#+ mod_validation, fig.width = 9, fig.height = 20, dpi = 300, out.width = '80%'
mod_val <- support_for_modules(mod_samples, modules_at_ages, palette = pal)

layout_val <- c(
  area(1,1,4,4),
  area(5,1,7,3),
  area(8,1,9,2)
  )
wrap_plots(rev(mod_val$plot), design = layout_val, guides = "keep")

mod_val$mean_support

  
#' ### Ancestral states at internal nodes of Colias phylogeny
#' We also wanted to do a tradition ancestral state reconstruction (ASR), calculating  interaction probabilities at internal nodes of the Colias tree. And now that we have defined modules for the extant network, we can also use them to group hosts in the ASR. 
#' 
#' **Interaction probability at internal nodes**
#' 

#+ ancestral_states, fig.width = 15, fig.height = 13, dpi = 300
at_nodes <- posterior_at_nodes(history, tree, host_tree)
p_asr <- plot_matrix_phylo(matrix, at_nodes, tree, host_tree, modules = mod, threshold = 0.9, colors = pal_extant)
p_asr[[2]] <- p_asr[[2]] + theme(axis.text = element_text(face = "italic"))
p_asr

#+ mrcas, echo = TRUE
# root state
root_state <- sort(at_nodes$post_states["Index_99",,], decreasing = T)
root_state <- root_state[which(root_state > 0.8)]
root_state

# Clade A
cladeA_state <- at_nodes$post_states[paste0("Index_", 51:56),,]
cladeA_state <- cladeA_state[,colSums(cladeA_state) > 0.8]
cladeA_state

mrcaA_state <- sort(at_nodes$post_states["Index_56",,], decreasing = T)
mrcaA_state <- mrcaA_state[which(mrcaA_state > 0.8)]
mrcaA_state

# Clade B MRCA
mrcaB_state <- sort(at_nodes$post_states["Index_98",,], decreasing = T)
mrcaB_state <- mrcaB_state[which(mrcaB_state > 0.8)]
mrcaB_state


#' ### Network structure
#' 

net <- replace(matrix, which(matrix[,] == 2), 1)
index_at_ages_summary(list("0" = net), index = "Q")
index_at_ages_summary(list("0" = net), index = "NODF")
visweb(net, type = "nested")


/*

#    N_obs      mean         sd         z age
# 35.90436  10.23888  0.6518014  39.37623   0
#  
#     Q_obs      mean        sd         z age
# 0.3998386 0.4311886 0.0113634 -2.758858   0

  
###
  
# Figures for paper

p_asr

p_mod_matrix_ages[[1]] +theme(axis.text = element_text(face = "italic"),
                              legend.position = "bottom",
                              title = element_blank())
  
plot(tree, show.node.label = TRUE, cex = 0.8)

plot_param

###  
*/







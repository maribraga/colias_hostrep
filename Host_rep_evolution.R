#'---
#'title: "Colias host repertoire evolution"
#'author: "Mariana P Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#+ setup, include = FALSE
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

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
log1 <- read.table(paste0(path_out,"out.2.2b.colias.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.3.2b.colias.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")

#' **Convergence test**
#' 
its <- seq(50000,500000,1000)

chain1 <- filter(log1, iteration %in% its) %>% as.mcmc()
chain2 <- filter(log2, iteration %in% its) %>% as.mcmc()

gelman.diag(mcmc.list(chain1, chain2))

#' All values are very close to 1, so we are good to go.
#' 

#' **Parameter estimates**
#' 
#' Now, we can plot the posterior distributions. Since both chains have samples from the same posterior distribution, we can use only one of the log files for that.
#' 
#+ densities, out.width = '50%', dpi = 300
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
history <- read_history(paste0(path_out, "out.2.2b.colias.history.txt"), burnin = 0) %>% 
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
mod <- readRDS(paste0(path_evol, "5mods.rds"))


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

pt_rev + geom_vline(xintercept = c(-2.1,-1.4,-0.7,0), col = "blue")

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

# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree, extant_modules = mod)

# modules and submodules
mods_pal <- sort(unique(modules_at_ages$matched_modules$original_and_matched_module_names$module_name))
pal <- c("#edae49", "#d1495b", "#d86357", "#df7c52", 
          "#9d5568", "#66a182", "#00798c", "#2e4057")
pal_extant <- c("#edae49", "#d1495b", "#66a182", "#00798c", "#2e4057")


#' **Plot ancestral networks**
#' 

#+ ancestral_nets, fig.width = 12, fig.height = 10, dpi = 300
p_nets <- plot_ancestral_networks(summary_networks, modules_at_ages, tree, palette = pal, node_size = 2)
wrap_plots(p_nets, heights = c(2,3), guides = "collect")



#' **Plot ancestral networks as matrices**
#' 
#' An alternative way to visualize the networks is by plotting them as matrices.
#' 

#+ anc_nets_matrix, fig.width = 7, fig.height = 22, dpi = 300
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

#+ mod_validation
samples_at_ages <- get_sampled_networks(at_ages)

#+ mod_val, eval = FALSE
mod_samples <- modules_from_samples(samples_at_ages)
saveRDS(mod_samples, "ignore/2s_new_tree/evolnets/mod_samples.rds")
# took 56 min!

#+ mod_val_read, echo = FALSE
mod_samples <- readRDS("ignore/2s_new_tree/evolnets/mod_samples.rds")

#+ support, fig.width = 7, fig.height = 14, dpi = 300
mod_val <- support_for_modules(mod_samples, modules_at_ages, palette = pal)
mod_val$mean_support

layout_val <- c(
  area(1,1,4,4),
  area(5,1,7,3),
  area(8,1,9,2)
  )
wrap_plots(rev(mod_val$plot), design = layout_val, guides = "collect")


  
#' ### Ancestral states at internal nodes of Colias phylogeny
#' We also wanted to do a tradition ancestral state reconstruction (ASR), calculating  interaction probabilities at internal nodes of the Colias tree. And now that we have defined modules for the extant network, we can also use them to group hosts in the ASR. 
#' 
#' **Interaction probability at internal nodes**
#' 

#+ ancestral_states, fig.width = 13, fig.height = 11, dpi = 300
at_nodes <- posterior_at_nodes(history, tree, host_tree)

p_asr <- plot_matrix_phylo(matrix, at_nodes, tree, host_tree, modules = match_mod_ext, threshold = 0.9, colors = pal_extant)
p_asr

#' **Plot ancestral states as matrix**
#' 

#+ nodes_matrix, fig.width = 12, fig.height = 6, dpi = 300
pp <- at_nodes$post_states[,,1]

graph <- igraph::graph_from_incidence_matrix(pp, weighted = TRUE)

nodes <- paste0("Index_",(Ntip(tree)+1):(Ntip(tree)+Nnode(tree)))

edge_list_nodes_mod <- igraph::get.data.frame(graph, what = "edges") %>% 
  dplyr::mutate(from = factor(from, levels = nodes),
                to = factor(to, levels = host_tree$tip.label),
                name = to) %>% 
  left_join(filter(match_mod_ext, type == "host")) %>% 
  rename(p = weight)

gg_all_nodes_mod <- ggplot(edge_list_nodes_mod, aes(x = to, y = from)) +
  geom_tile(aes(fill = module, alpha = p)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_discrete(type = pal_extant) +
  scale_alpha(guide = guide_legend(title = "Posterior\nprobability")) +
  labs(fill = "Module") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

gg_high_nodes_mods <- filter(edge_list_nodes_mod, p >= 0.9) %>%
  mutate(group = case_when(p < 0.95 ~ ".90-.94",
                           p >= 0.95 ~ ".95-1",
                           TRUE ~ as.character(p))) %>%
  ggplot(aes(x = to, y = from)) + 
  geom_tile(aes(fill = module, alpha = group)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_discrete(type = pal_extant[c(2,4,5)]) +
  scale_alpha_discrete(range = c(0.6,1), guide = guide_legend(title = "Posterior\nprobability")) +
  labs(fill = "Module") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

gg_all_nodes_mod + gg_high_nodes_mods

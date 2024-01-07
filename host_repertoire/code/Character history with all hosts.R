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


path_data <- "host_repertoire/data/"
path_evol <- "host_repertoire/evolnets/"
path_out <- "host_repertoire/output/"

matrix <- read.csv(paste0(path_data,"interaction_matrix_all_hosts.csv"), row.names = 1) %>% as.matrix()

#mod <- mycomputeModules(matrix)
#saveRDS(mod4, paste0(path_evol, "11mods.rds"))
mod <- readRDS(paste0(path_evol, "11mods.rds"))
plotModuleWeb(mod, labsize = 0.5)

tree <- read_tree_from_revbayes(paste0(path_evol,"tree_Rev.tre"))
host_tree <- phytools::starTree(colnames(matrix))
plot(host_tree)

hist5 <- read_history(paste0(path_out, "out.5_thin.txt"), burnin = 0) %>% 
   filter(iteration > 50000)
hist4 <- read_history(paste0(path_out, "out.4_thin.txt"), burnin = 0) %>% 
   filter(iteration > 50000)
hist2 <- read_history(paste0(path_out, "out.2_thin.txt"), burnin = 0) %>% 
   filter(iteration > 20000)

count_events(hist2)
count_events(hist4)
count_events(hist5)

at_nodes5 <- posterior_at_nodes(hist5, tree, host_tree)
at_nodes4 <- posterior_at_nodes(hist4, tree, host_tree)
at_nodes2 <- posterior_at_nodes(hist2, tree, host_tree)

pp5 <- at_nodes5$post_repertoires[,,2] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "node") %>% 
  pivot_longer(2:63, names_to = "host", values_to = "pp5")
  
pp4 <- at_nodes4$post_repertoires[,,2] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "node") %>% 
  pivot_longer(2:63, names_to = "host", values_to = "pp6")

pp2 <- at_nodes2$post_repertoires[,,2] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "node") %>% 
  pivot_longer(2:63, names_to = "host", values_to = "pp7")

pps <- pp5 %>% left_join(pp4) %>% left_join(pp2) %>%
  unite("pair", 1:2, remove = TRUE) %>% 
  pivot_longer(2:4, names_to = "chain", values_to = "pp")

ggplot(pps, aes(chain, pair)) +
  geom_raster(aes(alpha = pp))  +
  theme_bw() +
  theme(axis.text.y = element_blank())

p_asr5 <- plot_matrix_phylo(matrix, at_nodes5, tree, host_tree, modules = mod, threshold = 0.9)
p_asr5[[2]] <- p_asr5[[2]] + theme(axis.text = element_text(face = "italic", size = 6))
p_asr5[[4]] <- NULL
p_asr5  

p_asr4 <- plot_matrix_phylo(matrix, at_nodes4, tree, host_tree, modules = mod, threshold = 0.9)
p_asr4[[2]] <- p_asr4[[2]] + theme(axis.text = element_text(face = "italic", size = 6))
p_asr4[[4]] <- NULL
p_asr4  

p_asr2 <- plot_matrix_phylo(matrix, at_nodes2, tree, host_tree, modules = mod, threshold = 0.9, ladderize = FALSE)
p_asr2[[2]] <- p_asr2[[2]] + 
  scale_x_discrete(limits = host_tree$tip.label[mod@orderB]) +
  theme(axis.text = element_text(face = "italic", size = 6))
p_asr2[[4]] <- NULL
p_asr2 

p_asr2[[1]] + p_asr4[[1]]


#### Parameters ####

log1 <- read.table(paste0(path_out,"out.5.colias_m0_62h.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.7.colias_m0_62h.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:7)]
log2 <- log2[,c(1,5:7)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock", "gain", "loss")

its <- seq(4100,11000,100)

chain1 <- filter(log1, iteration %in% its)
chain2 <- filter(log2, iteration %in% its)

gelman.diag(mcmc.list(as.mcmc(chain1), as.mcmc(chain2)))

effectiveSize(chain1)
effectiveSize(chain2)
chain2 <- filter(log2, iteration > 1000)
effectiveSize(chain2)

## With beta

logbeta <- read.table(paste0(path_out,"out.6.colias.log"), header = TRUE)

logbeta <- logbeta[,c(1,5,7:9)]
colnames(logbeta) <- c("iteration","clock","beta", "gain", "loss")
chainbeta <- filter(logbeta, iteration > 5000)
effectiveSize(chainbeta)

## Estimates

parameters <- chain2 %>% 
  pivot_longer(cols = 2:4, names_to = "parameter", values_to = "value")

parameters %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))

ggplot(parameters) +
  geom_density(aes(value, fill = parameter, col = parameter, y = after_stat(scaled)), alpha = 0.6 ) +
  theme_bw()


param_beta <- chainbeta %>% 
  pivot_longer(cols = 2:5, names_to = "parameter", values_to = "value")

param_beta %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))


## Bayes factor - beta

d_prior <- dexp(x=0, rate=1)

kd_beta <- kdensity(x = chainbeta$beta, 
                    kernel='gamma', 
                    support=c(0,Inf), 
                    bw = 0.02)
max = kd_beta(0)

(BF <- d_prior/max)

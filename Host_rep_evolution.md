Colias host repertoire evolution
================
Mariana P Braga
05 August, 2022

------------------------------------------------------------------------

Script for analyses performed in … *…*

## Set up

For this script we’ll need *evolnets* to analyze the output posterior
distribution for character history. You can install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("maribraga/evolnets")
```

We also need other packages

``` r
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
```

## Parameter estimates

The first thing we need to do is to check that independent MCMC chains
have converged to the same posterior distribution. Both chains have
achieved an effective sample size ESS \> 200 for every parameter.

**log files**

``` r
# read files
path_out<- "ignore/2s_new_tree/output/"
log1 <- read.table(paste0(path_out,"out.5.2b.colias.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.2.2b.colias.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")
```

**Convergence test**

``` r
its <- seq(20000,200000,500)

chain1 <- filter(log1, iteration %in% its)
chain2 <- filter(log2, iteration %in% its)

gelman.diag(mcmc.list(as.mcmc(chain1), as.mcmc(chain2)))
```

    ## Potential scale reduction factors:
    ## 
    ##           Point est. Upper C.I.
    ## iteration        NaN        NaN
    ## clock          0.997      0.998
    ## beta           1.007      1.037
    ## gain           1.017      1.082
    ## loss           1.017      1.082
    ## 
    ## Multivariate psrf
    ## 
    ## 1.02

``` r
effectiveSize(chain1)
```

    ## iteration     clock      beta      gain      loss 
    ##    0.0000  278.8117  361.0000  361.0000  361.0000

All values are very close to 1, so we are good to go.

**Parameter estimates**

Now, we can plot the posterior distributions. Since both chains have
samples from the same posterior distribution, we can use only one of the
log files for that.

``` r
parameters <- chain1 %>% 
  pivot_longer(cols = 2:5, names_to = "parameter", values_to = "value")

ggplot(parameters, aes(parameter, value)) +
  geom_violin() +
  stat_summary(fun.data = "median_hilow", color = "red", size = 0.5) +
  theme_bw()
```

<img src="Host_rep_evolution_files/figure-gfm/densities-1.png" width="50%" />

``` r
parameters %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value),
            sd = sd(value))
```

    ## # A tibble: 4 × 3
    ##   parameter   mean      sd
    ##   <chr>      <dbl>   <dbl>
    ## 1 beta      0.0885 0.0834 
    ## 2 clock     0.586  0.0541 
    ## 3 gain      0.0694 0.00810
    ## 4 loss      0.931  0.00810

**Bayes factor**

The parameter called `beta` defines whether the phylogenetic distance
between hosts affects the probability of gaining new hosts. When
`beta = 0`, phylogenetic distances between hosts do not matter, hence
all hosts are equally likely to be gained. We use Bayes factor to test
whether the inferred `beta` is significantly different from 0. If the
factor is \< 1, it means that the phylogenetic distances do not matter
during host gain events.

``` r
d_prior <- dexp(x=0, rate=1)

kd_beta <- kdensity(x = log1$beta, 
                    kernel='gamma', 
                    support=c(0,Inf), 
                    bw = 0.02)
max = kd_beta(0)

(BF <- d_prior/max)
```

    ## [1] 0.1211106

## Character history

Let’s move on to reconstruction of the history of evolution of host
repertoire across the Colias phylogeny.

#### Data

**Trees**

We use `read_tree_from_revbayes` to read the Colias tree because this
file was exported from RevBayes and contains the node labels given by
RevBayes. This will be very important in the analysis!

``` r
path_data <- "ignore/2s_new_tree/data/"
path_evol <- "ignore/2s_new_tree/evolnets/"

tree <- read_tree_from_revbayes(paste0(path_evol,"tree_final_Rev.tre"))
host_tree <- read.tree(paste0(path_data,"host_tree.tre"))
```

**Extant network**

This matrix contains 0s and 2s because in the host repertoire model in
RevBayes, there are 3 possible states (0,1,2), where 1 means “potential
host” and 2 means “actual host”. We used the 2-state model for the
reconstruction in RevBayes, so we are only interested in the 0s and 2s,
no potential host.

``` r
matrix <- read.csv(paste0(path_data,"matrix_phylo_timetree.csv"), row.names = 1) %>% as.matrix()
```

**Read in .history.txt files**

We’ll use the *evolnets* function `read_history()` to read a file
outputed from RevBayes with sampled histories during MCMC

``` r
history <- read_history(paste0(path_out, "out.5.2b.colias.history.txt"), burnin = 0) %>% 
  filter(iteration %in% its)
```

#### Number of events and effective rate of evolution

``` r
# Estimated number of events across the Colias phylogeny
count_events(history)
```

    ##       mean HPD95.lower HPD95.upper
    ## 1 291.1025         251         345

``` r
# How many events were host gains and how many host losses?
gl <- count_gl(history)
# percentage of gains and losses
gl/sum(gl)
```

    ##     gains    losses 
    ## 0.4042326 0.5957674

``` r
# Considering the number of events and the total length of the Colias phylogeny, 
# what is the estimated rate of host repertoire evolution?
effective_rate(history, tree)
```

    ##       mean HPD95.lower HPD95.upper
    ## 1 4.186265    3.609562    4.961351

#### Modules of the extant (present-day) network

``` r
# find modules
mod <- mycomputeModules(matrix)
```

### Extant and ancestral networks

**Time points of interest (ages)**

The first step to reconstruct ancestral networks is to define the time
points during Colias diversification that we want to look at.

``` r
# visually determine interesting time points to reconstruct ancestral networks
pt <- ggtree(tree) + 
  geom_tiplab() + 
  geom_nodelab(size = 3, color = "grey40") + 
  theme_tree2()

pt_rev <- revts(pt) + xlim(c(-3, 0.8))

pt_rev + geom_vline(xintercept = c(-2.1,-1.4,-0.7,0), col = "blue")
```

<img src="Host_rep_evolution_files/figure-gfm/tree-1.png" width="60%" />

``` r
# choose time points
ages <- c(0,0.7,1.4,2.1)
```

I’ve chosen 2.1, 1.4, and 0.7, which means 2.1 Ma, 1.4 Ma and 700
thousand years ago.

**Interaction probability at given ages**

Now we calculate the posterior probability of interaction between each
host and each extant butterfly at each age in `ages`.

``` r
at_ages <- posterior_at_ages(history, ages, tree, host_tree)
```

**Summarize probabilities into ancestral networks**

Make binary or weighted networks? Discard interactions with posterior
probability \< threshold. We chose to reconstruct weighted networks with
a threshold of

``` r
summary_networks <- get_summary_networks(at_ages, threshold = 0.9)
```

``` r
# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree, extant_modules = mod)
```

``` r
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
```

**Plot ancestral networks**

``` r
p_nets <- plot_ancestral_networks(summary_networks, modules_at_ages, tree, palette = pal, node_size = 2)
wrap_plots(p_nets, heights = c(2,3), guides = "collect")
```

![](Host_rep_evolution_files/figure-gfm/ancestral_nets-1.png)<!-- -->

**Plot ancestral networks as matrices**

An alternative way to visualize the networks is by plotting them as
matrices.

``` r
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
```

<img src="Host_rep_evolution_files/figure-gfm/anc_nets_matrix-1.png" width="80%" />

**Module validation**

We use modules to facilitate visualisation but sometimes they also have
biological meaning. This step tells us whether the modules in the
ancestral networks are robust across MCMC samples or an artifact of the
summary networks.

``` r
samples_at_ages <- get_sampled_networks(at_ages)
```

``` r
mod_samples <- modules_from_samples(samples_at_ages)
```

``` r
mod_val <- support_for_modules(mod_samples, modules_at_ages, palette = pal)

layout_val <- c(
  area(1,1,4,4),
  area(5,1,7,3),
  area(8,1,9,2)
  )
wrap_plots(rev(mod_val$plot), design = layout_val, guides = "keep")
```

<img src="Host_rep_evolution_files/figure-gfm/mod_validation-1.png" width="80%" />

``` r
mod_val$mean_support
```

    ## $`2.1`
    ##   module      mean  geo_mean
    ## 1     M5 0.7377655 0.7034604
    ## 2     M6 0.6391136 0.5839256
    ## 
    ## $`1.4`
    ##   module      mean  geo_mean
    ## 1     M3 0.6246255 0.5563558
    ## 2     M4 0.8752597 0.8682388
    ## 3     M5 0.3634623 0.2377797
    ## 4     M6 0.5188920 0.3744296
    ## 
    ## $`0.7`
    ##   module      mean  geo_mean
    ## 1     M2 0.4938614 0.1676112
    ## 2   M3.1 0.5004693 0.4183086
    ## 3   M3.2 0.4972053 0.4042709
    ## 4     M4 0.7765755 0.7613942
    ## 5     M5 0.3885626 0.2253195
    ## 6     M6 0.3821607 0.2727022

### Ancestral states at internal nodes of Colias phylogeny

We also wanted to do a tradition ancestral state reconstruction (ASR),
calculating interaction probabilities at internal nodes of the Colias
tree. And now that we have defined modules for the extant network, we
can also use them to group hosts in the ASR.

**Interaction probability at internal nodes**

``` r
at_nodes <- posterior_at_nodes(history, tree, host_tree)
p_asr <- plot_matrix_phylo(matrix, at_nodes, tree, host_tree, modules = mod, threshold = 0.9, colors = pal_extant)
p_asr
```

![](Host_rep_evolution_files/figure-gfm/ancestral_states-1.png)<!-- -->

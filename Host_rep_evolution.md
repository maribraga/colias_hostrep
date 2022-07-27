Colias host repertoire evolution
================
Mariana P Braga
27 July, 2022

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
log1 <- read.table(paste0(path_out,"out.2.2b.colias.log"), header = TRUE)
log2 <- read.table(paste0(path_out,"out.3.2b.colias.log"), header = TRUE)

# take only columns of interest
log1 <- log1[,c(1,5:6,8:9)]
log2 <- log2[,c(1,5:6,8:9)]

# give them better column names
colnames(log1) <- colnames(log2) <- c("iteration","clock","beta", "gain", "loss")
```

**Convergence test**

``` r
its <- seq(50000,500000,1000)

chain1 <- filter(log1, iteration %in% its) %>% as.mcmc()
chain2 <- filter(log2, iteration %in% its) %>% as.mcmc()

gelman.diag(mcmc.list(chain1, chain2))
```

    ## Potential scale reduction factors:
    ## 
    ##           Point est. Upper C.I.
    ## iteration        NaN        NaN
    ## clock           1.01       1.03
    ## beta            1.00       1.00
    ## gain            1.00       1.02
    ## loss            1.00       1.02
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

All values are very close to 1, so we are good to go.

**Parameter estimates**

Now, we can plot the posterior distributions. Since both chains have
samples from the same posterior distribution, we can use only one of the
log files for that.

``` r
parameters <- log1 %>% 
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
    ## 1 beta      0.0964 0.0965 
    ## 2 clock     0.723  0.0648 
    ## 3 gain      0.0826 0.00907
    ## 4 loss      0.917  0.00907

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

    ## [1] 0.1247874

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
history <- read_history(paste0(path_out, "out.2.2b.colias.history.txt"), burnin = 0) %>% 
  filter(iteration %in% its)
```

#### Number of events and effective rate of evolution

``` r
# Estimated number of events across the Colias phylogeny
count_events(history)
```

    ##       mean HPD95.lower HPD95.upper
    ## 1 400.7938         337         462

``` r
# How many events were host gains and how many host losses?
gl <- count_gl(history)
# percentage of gains and losses
gl/sum(gl)
```

    ##     gains    losses 
    ## 0.4235553 0.5764447

``` r
# Considering the number of events and the total length of the Colias phylogeny, 
# what is the estimated rate of host repertoire evolution?
effective_rate(history, tree)
```

    ##       mean HPD95.lower HPD95.upper
    ## 1 5.763706    4.846305    6.643895

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

# find modules in extant and ancestral networks
modules_at_ages <- modules_across_ages(summary_networks, tree, extant_modules = mod)

# save names of matched modules for extant network for other plots
match_mod_ext <- modules_at_ages$matched_modules$nodes_and_modules_per_age %>% 
  filter(age == 0) %>% 
  select(name, module, type)

# modules and submodules
mods_pal <- sort(unique(modules_at_ages$matched_modules$original_and_matched_module_names$module_name))
pal <- c("#edae49", "#d1495b", "#d86357", "#df7c52", 
          "#9d5568", "#66a182", "#00798c", "#2e4057")
pal_extant <- c("#edae49", "#d1495b", "#66a182", "#00798c", "#2e4057")
```

``` r
mod_order <- match_mod_ext %>% filter(type == "host") %>% arrange(module) %>% pull(name)

mod_order_para <- match_mod_ext %>% filter(type == "symbiont") %>% arrange(module) %>% pull(name)

p_ext_mods <- plot_extant_matrix(matrix, match_mod_ext, 
                   parasite_order = mod_order_para, 
                   host_order = mod_order)
p_ext_mods + scale_fill_manual(values = pal_extant)
```

<img src="Host_rep_evolution_files/figure-gfm/mod_ext_matrix-1.png" width="70%" />

**Plot networks**

``` r
p_nets <- plot_ancestral_networks(summary_networks, modules_at_ages, tree, palette = pal, node_size = 2)
wrap_plots(p_nets, heights = c(2,3), guides = "collect")
```

![](Host_rep_evolution_files/figure-gfm/ancestral_nets-1.png)<!-- -->

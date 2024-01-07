library(RevGadgets)
library(tidyverse)
library(ape)

tree_file = "./biogeography/server_ignore/output_50k/output2/ase.tre"
log_file = "./biogeography/server_ignore/output_50k/output2/model.log"
states_file = "./biogeography/server_ignore/output_50k/output2/states.png"
posterior_file = "./biogeography/server_ignore/output_50k/output2/posterior.png"

## find out code given to states ##

tree <- read.tree("./biogeography/data/tree_colias_geo.tre")
tips <- tree$tip.label
plot(tree, cex = 0.5)

original_states <- read.csv("./biogeography/data/colias_geo_tbl.csv") %>% 
  filter(species %in% tips) %>% 
  arrange(factor(species, levels = tips))

tree$tip.label <- pull(original_states, range)
plot(tree, cex = 0.5)

labels <- c("1" = "EP", "2" = "WP", "3" = "NA", "4" = "PTP", "5" = "AF", "6" = "NT", 
            "9" = "PaLo", "11" = "NaA", "15" = "Wptp", "31" = "Arc", "32" = "EAs", "34" = "Pal",
            "14" = "Eptp", "46" = "PalAf", "67" = "North", "108" = "All", "25" = "Americas")

colors <- c("#00A2FF","#013459","#DA3541","#F8C700","#EF9FBA","#FF8C01",
            "#347FD5","#AA74A1","#097C65","#6E4672","#00D164","#009F52",
            "grey20","grey40","grey60","grey80", "darkred")
names(colors) <- labels

# from error message - vary across runs
all_states <- c(0, 1, 108, 11, 12, 13, 14, 15, 16, 19, 2, 20, 21, 23, 24, 25, 26, 27, 3, 31, 32, 34, 4, 46, 5, 56, 6, 67, 74, 84, 87, 9)

all_states <- as.character(sort(all_states))
extra_states <- setdiff(all_states, names(labels))
names(extra_states) <- extra_states

labels_all <- c(labels, extra_states)

states <- processAncStates(tree_file, state_labels = labels_all)

plotAncStatesMAP(t=states,
                 tree_layout="circular",
                 tip_labels_offset = 0.4,
                 node_size=1.5,
                 node_color_as="state",
                 node_color=colors,
                 node_size_as=NULL) +
                 ggplot2::theme(legend.position="bottom",
                                legend.title=element_blank())

ggsave(states_file, width = 9, height = 9)

trace <- readTrace(path=log_file)
plotTrace(trace, vars=c("phi_d[1]"))

ggsave(posterior_file, width = 9, height = 9)

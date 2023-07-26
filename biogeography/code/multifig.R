library(RevGadgets)
library(tidyverse)
library(ape)

tree_file = "./biogeography/server/output1/ase.tre"
log_file = "./biogeography/server/output1/model.log"
states_file = "./biogeography/server/output1/states.png"
posterior_file = "./biogeography/server/output1/posterior.png"

## find out code given to states ##

tree <- read.tree("./biogeography/data/tree_colias_geo.tre")
tips <- tree$tip.label
plot(tree, cex = 0.5)

original_states <- read.csv("./biogeography/data/Colias_distribution.csv") %>% 
  filter(species %in% tips) %>% 
  arrange(factor(species, levels = tips))

tree$tip.label <- pull(original_states, range)
plot(tree, cex = 0.5)

labels <- c("1" = "EP", "2" = "WP", "3" = "NA", "4" = "PTP", "5" = "AF", "6" = "NT", 
            "9" = "PaLo", "11" = "NaA", "15" = "Wptp", "31" = "Arc", "32" = "EAs", "34" = "Pal",
            "14" = "Eptp", "46" = "PalAf", "67" = "North", "108" = "All")

colors <- c("#00A2FF","#013459","#DA3541","#F8C700","#EF9FBA","#FF8C01",
            "#347FD5","#AA74A1","#097C65","#6E4672","#00D164","#009F52",
            "grey20","grey40","grey60","grey80")
names(colors) <- labels

# from error message
#all_states <- c(0, 1, 108, 11, 12, 13, 14, 15, 19, 2, 20, 21, 22, 25, 26, 27, 3, 31, 32, 34, 4, 46, 5, 52, 56, 6, 67, 74, 81, 86, 9, 96) # 1000gen
#all_states <- c(0, 1, 10, 108, 11, 12, 13, 14, 15, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 29, 3, 31, 32, 33, 34, 4, 40, 45, 46, 5, 52, 56, 6, 67, 74, 81, 86, 9, 96) # 20k out1
all_states <- c(0, 1, 108, 11, 12, 125, 13, 14, 15, 16, 19, 2, 20, 21, 23, 24, 25, 26, 27, 29, 3, 31, 32, 33, 34, 36, 4, 40, 42, 46, 5, 55, 56, 6, 67, 74, 81, 9) # 20k out2
all_states <- as.character(sort(all_states))
extra_states <- setdiff(all_states, names(labels))
names(extra_states) <- extra_states

labels_all <- c(labels, extra_states)

# labels <- c("1" = "EP", "2" = "WP", "3" = "NA", "4" = "PTP", "5" = "AF", "6" = "NT", "7" = "",
#             "8" = "", 9"" = "PaL", "10" = "", "11" = "NaA", "12" = "", "13" = "", "14" = "",
#             "15" = "Wptp", "16" = "", "17" = "", "18" = "", "19" = "", "20" = "", "21" = "",
#             "22" = "", "23" = "", "24" = "", "25" = "", "26" = "", "27" = "", "28" = "",
#             "29" = "", "30" = "", "31" = "Arc", "32" = "EAs", "33" = "", "34" = "Pal")

##

#labels <- c("0" = "Aa","1" = "Ca","2" = "Pa","3" = "Cc","4" = "Ad","5" = "El","6" = "AaCa","7" = "AaPa","8" = "CaPa","9" = "AaCc","10" = "CaCc","11" = "PaCc","12" = "AaAd","13" = "CaAd","14" = "PaAd","15" = "CcAd","16" = "AaEl","17" = "CaEl","18" = "PaEl","19" = "CcEl","20" = "AdEl","21" = "AaCaPa","22" = "AaCaCc","23" = "AaPaCc","24" = "CaPaCc","25" = "AaCaAd","26" = "AaPaAd","27" = "CaPaAd","28" = "AaCcAd","29" = "CaCcAd","30" = "PaCcAd","31" = "AaCaEl","32" = "AaPaEl","33" = "CaPaEl","34" = "AaCcEl","35" = "CaCcEl","36" = "PaCcEl","37" = "AaAdEl","38" = "CaAdEl","39" = "PaAdEl","40" = "CcAdEl","41" = "AaCaPaCc","42" = "AaCaPaAd","43" = "AaCaCcAd","44" = "AaPaCcAd","45" = "CaPaCcAd","46" = "AaCaPaEl","47" = "AaCaCcEl","48" = "AaPaCcEl","49" = "CaPaCcEl","50" = "AaCaAdEl","51" = "AaPaAdEl","52" = "CaPaAdEl","53" = "AaCcAdEl","54" = "CaCcAdEl","55" = "PaCcAdEl","56" = "AaCaPaCcAd","57" = "AaCaPaCcEl","58" = "AaCaPaAdEl","59" = "AaCaCcAdEl","60" = "AaPaCcAdEl","61" = "CaPaCcAdEl","62" = "AaCaPaCcAdEl")
#colors <- c("#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", "#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72")

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
library(RevGadgets)
library(tidyverse)
library(ape)
library(patchwork)

tree_file = "./biogeography/with_zerene/server_ignore/output/output2/ase.tre"
log_file = "./biogeography/with_zerene/server_ignore/output/output2/model.log"
states_file = "./biogeography/with_zerene/server_ignore/output/output2/states.pdf"
posterior_file = "./biogeography/with_zerene/server_ignore/output/output2/posterior.pdf"

## find out code given to states ##

tree <- read.tree("./biogeography/with_zerene/data/tree_colias_geo.tre")
tips <- tree$tip.label
plot(tree, cex = 0.5)

original_states <- read.csv("./biogeography/with_zerene/data/colias_geo_tbl.csv") %>% 
  arrange(factor(species, levels = tips))

tree$tip.label <- pull(original_states, range)
plot(tree, cex = 0.5)

labels <- c("1" = "EP", "2" = "WP", "3" = "NA", "4" = "PTP", "5" = "AF", "6" = "NT", 
            "9" = "PaLo", "11" = "NaA", "15" = "Wptp", "31" = "Arc", "32" = "EAs", "34" = "Pal",
            "14" = "Eptp", "46" = "PalAf", "67" = "North", "108" = "All", "25" = "Americas")

palette <- c("#00A2FF","#013459","#DA3541","#F8C700","#EF9FBA","#FF8C01",
            "#347FD5","#AA74A1","#097C65","#6E4672","#00D164","#009F52",
            "grey20","grey40","grey60","grey80", "darkred")
names(palette) <- labels

# from error message - vary across runs
all_states <- c(0, 1, 108, 11, 12, 13, 14, 15, 16, 19, 2, 20, 21, 23, 24, 25, 26, 27, 3, 31, 32, 34, 4, 46, 5, 56, 6, 67, 74, 84, 87, 9)

all_states <- as.character(sort(all_states))
extra_states <- setdiff(all_states, names(labels))
names(extra_states) <- extra_states

labels_all <- c(labels, extra_states)

states <- processAncStates(tree_file, state_labels = labels_all)

plotAncStatesMAP(t=states,
                 cladogenetic = TRUE,
                 tree_layout="rectangular",
                 tip_labels_offset = 0.1,
                 node_color_as="state",
                 node_color=palette,
                 state_transparency = 0.85,
                 #node_size = 2,
                 node_size_as="state_posterior",
                 timeline = TRUE,
                 geo_units = "epochs")

#ggsave(states_file, width = 9, height = 9)

## parameters
trace <- readTrace(path=log_file)
View(trace[[1]])

sigmas <- c("sigma_b[1]","sigma_b[2]","sigma_b[3]","sigma_d[1]","sigma_d[2]","sigma_d[3]","sigma_e[1]","sigma_e[2]","sigma_e[3]","sigma_w[1]","sigma_w[2]","sigma_w[3]")
phis <- c("phi_b[1]","phi_b[2]","phi_b[3]","phi_d[1]","phi_d[2]","phi_d[3]","phi_e[1]","phi_e[2]","phi_e[3]","phi_w[1]","phi_w[2]","phi_w[3]")
rhos <- c("rho_b","rho_d","rho_e","rho_w")

sigma_list <- list()
for(i in seq_along(sigmas)){
  sigma_list[i] <- plotTrace(trace, vars = sigmas[i])  
}
wrap_plots(sigma_list, ncol = 2)


phi_list <- list()
for(i in seq_along(phis)){
  phi_list[i] <- plotTrace(trace, vars = phis[i])  
}
wrap_plots(phi_list, ncol = 2)


rho_list <- list()
for(i in seq_along(rhos)){
  rho_list[i] <- plotTrace(trace, vars = rhos[i])  
}
wrap_plots(rho_list, ncol = 2)


plotTrace(trace, vars=c("phi_d[1]"))
#ggsave(posterior_file, width = 9, height = 9)

## speciation and extinction
bd <- processSSE(path=log_file, speciation = "lambda", extinction = "mu")
View(bd)

bd_subset <- filter(bd, observed_state %in% 0:6)
plotMuSSE(bd_subset)
bd_subset <- filter(bd, observed_state %in% 6:17)
plotMuSSE(bd_subset)


# joint posterior distribution

# 2d density

proc_groups <- c("_b", "_d", "_e", "_w")
phis <- array(data = paste0("phi",proc_groups), dim = c(4,1)) %>% 
  apply(1, FUN = function(x) paste0(x,"[",1:2,"]")) %>% 
  as.vector()
sigmas <- array(data = paste0("sigma",proc_groups), dim = c(4,1)) %>% 
  apply(1, FUN = function(x) paste0(x,"[",1:2,"]")) %>% 
  as.vector()

plot_2d <- function(trace, x = NULL, y = NULL){
  p <- ggplot(trace[[1]], aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(shape = 21, alpha = 0.5, colour = "white", fill = "blue", size = 1) +
    geom_hline(yintercept = 0, col = "red") +
    geom_vline(xintercept = 0, col = "red") +
    theme_bw()
  p
}

list_2d <- list()
for(i in seq_along(sigmas)){
  list_2d[[i]] <- plot_2d(trace, x = phis[i], y = sigmas[i])
}

wrap_plots(list_2d, ncol=2)
plot_2d_file <- here("biogeography/with_zerene/server_ignore/output/output2/2d_plots.pdf")
#ggsave(plot_2d_file, width = 9, height = 12)

# calculate joint pp

trace_d <- trace[[1]] %>% 
  select(Iteration, `phi_d[1]`, `sigma_d[1]`) %>% 
  mutate(joint_p = case_when(`phi_d[1]` >= 0 & `sigma_d[1]` >= 0 ~ "both_positive",
                             TRUE ~ "one_or_both_negative"))
trace_e <- trace[[1]] %>% 
  select(Iteration, `phi_e[1]`, `sigma_e[1]`) %>% 
  mutate(joint_p = case_when(`phi_e[1]` >= 0 & `sigma_e[1]` >= 0 ~ "both_positive",
                             TRUE ~ "one_or_both_negative"))
trace_b <- trace[[1]] %>% 
  select(Iteration, `phi_b[1]`, `sigma_b[1]`) %>% 
  mutate(joint_p = case_when(`phi_b[1]` >= 0 & `sigma_b[1]` >= 0 ~ "both_positive",
                             TRUE ~ "one_or_both_negative"))
trace_w <- trace[[1]] %>% 
  select(Iteration, `phi_w[1]`, `sigma_w[1]`) %>% 
  mutate(joint_p = case_when(`phi_w[1]` >= 0 & `sigma_w[1]` >= 0 ~ "both_positive",
                             TRUE ~ "one_or_both_negative"))

trace_d %>% 
  group_by(joint_p) %>% 
  summarise(n = n()/nrow(.))

trace_e %>% 
  group_by(joint_p) %>% 
  summarise(n = n()/nrow(.))

trace_b %>% 
  group_by(joint_p) %>% 
  summarise(n = n()/nrow(.))

trace_w %>% 
  group_by(joint_p) %>% 
  summarise(n = n()/nrow(.))


###### Overall effect of parameters over areas ####



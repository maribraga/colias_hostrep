library(RevGadgets)
library(tidyverse)
library(ape)
library(patchwork)
#library(MCMCpack)
library(coda)
library(here)

#### Analyses with mean altitude

log_file <- here("biogeography/without_zerene/server_ignore/output/run1/model.log")
log_file2 <- here("biogeography/without_zerene/server_ignore/output/run2/model.log")
geo_log_file <- here("biogeography/without_zerene/server_ignore/output/run2/geo_features.log")

tree_file <- here("biogeography/without_zerene/server_ignore/output/run2/ase.tre")
state_label_file <- here("biogeography/without_zerene/server_ignore/output/run2/state_labels.txt")
anc_states_file <- here("biogeography/without_zerene/server_ignore/output/run2/anc_states.pdf")
plot_2d_file <- here("biogeography/without_zerene/server_ignore/output/run2/2d_plots.pdf")

#### find out code given to states ####

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



#### speciation and extinction ####

bd <- processSSE(path=log_file2, speciation = "lambda", extinction = "mu")
View(bd)

bd_subset <- filter(bd, observed_state %in% 0:6)
plotMuSSE(bd_subset)
bd_subset <- filter(bd, observed_state %in% 6:17)
plotMuSSE(bd_subset)




###### Overall effect of parameters over areas ####

library(bayestestR)
library(ggraph)
library(tidygraph)

trace_geo <- readTrace(path = geo_log_file)
View(trace_geo[[1]])
# correct colnames
#colnames(trace_geo[[1]]) <- c("Iteration","Posterior","Likelihood","Prior","m_b[1][1]","m_b[1][2]","m_b[1][3]","m_b[1][4]","m_b[1][5]","m_b[1][6]","m_b[1][7]","m_b[2][1]","m_b[2][2]","m_b[2][3]","m_b[2][4]","m_b[2][5]","m_b[2][6]","m_b[2][7]","m_b[3][1]","m_b[3][2]","m_b[3][3]","m_b[3][4]","m_b[3][5]","m_b[3][6]","m_b[3][7]","m_b[4][1]","m_b[4][2]","m_b[4][3]","m_b[4][4]","m_b[4][5]","m_b[4][6]","m_b[4][7]","m_b[5][1]","m_b[5][2]","m_b[5][3]","m_b[5][4]","m_b[5][5]","m_b[5][6]","m_b[5][7]","m_b[6][1]","m_b[6][2]","m_b[6][3]","m_b[6][4]","m_b[6][5]","m_b[6][6]","m_b[6][7]","m_b[7][1]","m_b[7][2]","m_b[7][3]","m_b[7][4]","m_b[7][5]","m_b[7][6]","m_b[7][7]","m_d[1][1]","m_d[1][2]","m_d[1][3]","m_d[1][4]","m_d[1][5]","m_d[1][6]","m_d[1][7]","m_d[2][1]","m_d[2][2]","m_d[2][3]","m_d[2][4]","m_d[2][5]","m_d[2][6]","m_d[2][7]","m_d[3][1]","m_d[3][2]","m_d[3][3]","m_d[3][4]","m_d[3][5]","m_d[3][6]","m_d[3][7]","m_d[4][1]","m_d[4][2]","m_d[4][3]","m_d[4][4]","m_d[4][5]","m_d[4][6]","m_d[4][7]","m_d[5][1]","m_d[5][2]","m_d[5][3]","m_d[5][4]","m_d[5][5]","m_d[5][6]","m_d[5][7]","m_d[6][1]","m_d[6][2]","m_d[6][3]","m_d[6][4]","m_d[6][5]","m_d[6][6]","m_d[6][7]","m_d[7][1]","m_d[7][2]","m_d[7][3]","m_d[7][4]","m_d[7][5]","m_d[7][6]","m_d[7][7]","m_e[1]","m_e[2]","m_e[3]","m_e[4]","m_e[5]","m_e[6]","m_e[7]","m_w[1]","m_w[2]","m_w[3]","m_w[4]","m_w[5]","m_w[6]","m_w[7]","rho_b","rho_d","rho_e","rho_w")



# between-region speciation

mbs <- grep("m_b", colnames(trace_geo[[1]]), value = TRUE)
keep <- str_extract_all(mbs, "\\d") %>% lapply(FUN = function(x) x[1] != x[2]) %>% unlist()
mbs <- mbs[keep]
mbs 

mb_maps <- as_tibble(map_estimate(trace_geo[[1]][,mbs])) %>% 
  rename(mb = Parameter, weight = MAP_Estimate) %>% 
  mutate(from = str_extract(mb, "\\d") %>% as.integer(),
         to = str_extract_all(mb, "\\d", simplify = TRUE)[,2]  %>% as.integer())

mb_tbl <- mb_maps %>% 
  select(from, to, weight)


# within-region speciation

mws <- grep("m_w", colnames(trace_geo[[1]]), value = TRUE)
mw_tbl <- as_tibble(map_estimate(trace_geo[[1]][,mws])) %>% 
  rename(area = Parameter, mw = MAP_Estimate) %>%
  mutate(index = 1:length(mws),
         name = state_labels$combination[1:length(mws)]) %>% 
  select(index, name, mw)

graph <- tbl_graph(nodes = mw_tbl, edges = mb_tbl)
graph


palette_fill <- c("NTr" = "#FF8C01",
                  "NAr" = "#DA3541",
                  "EPa" = "#00A2FF",
                  "PTH" = "#F8C700",
                  "WPa" = "#013459",
                  "Afr" = "#EF9FBA",
                  "SEA" = "white")

palette_col <- c("NTr" = "#FF8C01",
                 "NAr" = "#DA3541",
                 "EPa" = "#00A2FF",
                 "PTH" = "#F8C700",
                 "WPa" = "#013459",
                 "Afr" = "#EF9FBA",
                 "SEA" = "grey30")

ggraph(graph, layout = 'stress') +  
  geom_edge_link(aes(edge_width = weight)) + 
  geom_node_point(aes(color = name, size = mw)) +
  scale_color_manual(values = palette_col) +
  theme_void()



# between-region dispersal

mds <- grep("m_d", colnames(trace_geo[[1]]), value = TRUE)
keep <- str_extract_all(mds, "\\d") %>% lapply(FUN = function(x) x[1] != x[2]) %>% unlist()
mds <- mds[keep]
mds 

md_maps <- as_tibble(map_estimate(trace_geo[[1]][,mds])) %>% 
  rename(md = Parameter, weight = MAP_Estimate) %>% 
  mutate(from = str_extract(md, "\\d") %>% as.integer(),
         to = str_extract_all(md, "\\d", simplify = TRUE)[,2]  %>% as.integer())

md_tbl <- md_maps %>% 
  select(from, to, weight)


# within-region extinction

mes <- grep("m_e", colnames(trace_geo[[1]]), value = TRUE)
me_tbl <- as_tibble(map_estimate(trace_geo[[1]][,mes])) %>% 
  rename(area = Parameter, me = MAP_Estimate) %>%
  mutate(index = 1:length(mes),
         name = state_labels$combination[1:length(mes)]) %>% 
  select(index, name, me)

graph2 <- tbl_graph(nodes = me_tbl, edges = mb_tbl, directed = TRUE)
graph2

ggraph(graph2, layout = 'stress') +  
  geom_edge_link(aes(edge_width = weight)) + 
  geom_node_point(aes(color = name, size = me)) +
  scale_color_manual(values = palette_col) +
  theme_void()


# Sum of rates out and in of each region

m_de <- md_tbl %>% 
  left_join(select(me_tbl, -me), by = c("from" = "index")) %>% 
  rename(from_name = name) %>% 
  left_join(select(me_tbl, -me), by = c("to" = "index")) %>% 
  rename(to_name = name, dispersal = weight)

sum_out_mde <- m_de %>% 
  group_by(from_name) %>% 
  summarise(sum_out = sum(dispersal)) %>% 
  arrange(sum_out)

sum_md <- m_de %>% 
  group_by(to_name) %>% 
  summarise(sum_in = sum(dispersal)) %>% 
  arrange(sum_in) %>% 
  left_join(sum_out_mde, by = c("to_name" = "from_name"))

sum_md

pmd <- ggplot(m_de) +
  geom_raster(aes(from_name, to_name, alpha = dispersal)) +
  theme_bw()

pme <- ggplot(me_tbl) +
  geom_col(aes(name, me, fill = name)) +
  scale_fill_manual(values = palette_col) +
  theme_bw()

pmd_in <- ggplot(sum_md) +
  geom_col(aes(to_name, sum_in, fill = to_name)) +
  scale_fill_manual(values = palette_col) +
  theme_bw()

pmd_out <- ggplot(sum_md) +
  geom_col(aes(to_name, sum_out, fill = to_name)) +
  scale_fill_manual(values = palette_col) +
  theme_bw()


# Sum of speciation rates for each region

m_bw <- mb_tbl %>% 
  left_join(select(mw_tbl, -mw), by = c("from" = "index")) %>% 
  rename(from_name = name) %>% 
  left_join(select(mw_tbl, -mw), by = c("to" = "index")) %>% 
  rename(to_name = name, between = weight)

sum_out_mb <- m_bw %>% 
  group_by(from_name) %>% 
  summarise(sum_out = sum(between)) %>% 
  arrange(sum_out)

sum_mb <- m_bw %>% 
  group_by(to_name) %>% 
  summarise(sum_in = sum(between)) %>% 
  arrange(sum_in) %>% 
  left_join(sum_out_mb, by = c("to_name" = "from_name"))

sum_mb

pmb <- ggplot(m_bw) +
  geom_raster(aes(from_name, to_name, alpha = between)) +
  theme_bw()

pmw <- ggplot(mw_tbl) +
  geom_col(aes(name, mw, fill = name)) +
  scale_fill_manual(values = palette_col) +
  theme_bw()

(pmd + pme) / (pmb + pmw)

(pmd_in + pmd_out) / (pme + pmw)


#### Altitude variance instead of mean

log_file <- here("biogeography/without_zerene/server_ignore/output/altvar/run1/model.log")
log_file2 <- here("biogeography/without_zerene/server_ignore/output/altvar/run2/model.log")
tree_file <- here("biogeography/without_zerene/server_ignore/output/altvar/run2/ase.tre")
geo_log_file <- here("biogeography/without_zerene/server_ignore/output/altvar/run2/geo_features.log")


# install branch from RPanda repo
#install_github("hmorlon/PANDA", ref = "clade.shift.model")

# load packages
library(ape)
library(RPANDA)
library(RColorBrewer)


setwd("~/repos/colias_hostrep/diversification")

tree <- readRDS("tree_Colias.rds") 
groups <- readRDS("groups.rds")

sfrac <- get.sampling.fractions(tree, groups, plot = TRUE, cex = 0.3)
node_comb <- get.comb.shift(tree, groups, sfrac)

shifts <- shift.estimates(phylo = tree,
                          data = groups,
                          sampling.fractions = sfrac,
                          comb.shift = node_comb,
                          models = c("BCST","BVAR",
                                     "BCST_DCST","BVAR_DVAR",
                                     "BVAR_DCST","BCST_DVAR"),
                          backbone.option = "crown.shift",
                          multi.backbone = "all",
                          Ncores = 4)

head(shifts$total)
shifts$backbones$`75/`
shifts$subclades$`75`

plot.phylo.comb(phylo = tree, 
                data = groups, 
                sampling.fractions = sfrac,
                shift.res = shifts,
                combi = 1, 
                edge.width = 2, 
                label.offset = 0.3, 
                lad = F, 
                cex = 0.4)
axisPhylo()

rates <- div.rates(phylo = tree,
                   shift.res = shifts,
                   time.interval = 0.5,
                   combi = 1, part = "all")
rates

paleodiversity <- paleodiv(phylo = tree,
                           data = groups,
                           sampling.fractions = sfrac,
                           shift.res = shifts,
                           time.interval = 0.5,
                           split.div = TRUE)
paleodiversity

crown_age <- max(node.depth.edgelength(tree))
time <- -c(crown_age, seq(floor(crown_age), 0, by = -0.5))
plot(time, paleodiversity["backbone",],
     bty = "n", las = 1, ylim = c(0, 80),
     main = "Paleodiversity dynamic of Colias",
     xlab = "Time (Myrs)", ylab = "Number of species", lwd = 2, type = "l", lty = 1)


#### PLOT everything together ----

# Define the layout of plot with the matrix mat 
mat <- matrix(c(1,1,2,3), nrow = 2, ncol = 2)

layout(mat) 
par(mar = c(6,4,3,4)) 

subclade_colors <- brewer.pal(8, "Dark2")[1:4] 
group_colors <- c(subclade_colors, "black")

# Phylogeny with shifts 

plot.phylo.comb(phylo = tree, data = groups, sampling.fractions = sfrac, shift.res = shifts, combi = 1, label.offset = 0.1, cex = 0.7, lad = F) 
axisPhylo() 
mtext("Time (Ma)", side = 1, line = 2, cex = 1, at = 17)

# Rates through time 
#par(mar = c(6,4,3,4)) 
plot(time, rates[[1]][1,], type = "l", col = group_colors[1], las = 1, xlab = "Time (Myrs)", ylab = "Rates (Events/ Lineages / Myr)", main = "Diversification rates")
lines(time, rates[[length(rates)]][1,], type = "l", col = "black") 
legend("topright", legend = c("Backbone", "After shift"), col = c("black", group_colors[1]), lty = 1)

# Paleodiversity dynamics 
#par(mar = c(6,4,3,4)) 
plot(time, paleodiversity["75",],
     las = 1, ylim = c(0, 80),
     col = group_colors[1],
     main = "Paleodiversity dynamics of Colias",
     xlab = "Time (Myrs)", ylab = "Number of species", type = "l") 
lines(time, paleodiversity["backbone",], col = "black")
legend("topleft", legend = c("Backbone", "After shift"), col = c("black", group_colors[1]), lty = 1)


#### Model adequacy ----

all_posteriors <- simul.comb.shift(phylo = tree, sampling.fractions = sfrac, shift.res = shifts) 
posterior_trees <- all_posteriors #[1:500]

# Number of tips per group

ntip_by_groups <- lapply(posterior_trees, function(x) table(sapply(strsplit(x$tip.label, split = ""), "[[", 1))) 
ntip_by_groups_df <- as.data.frame(do.call(rbind, ntip_by_groups)) 
ntip_by_groups_colias <- c(73,7)


par(mfrow = c(1,2), cex = 0.6) 
boxplot(ntip_by_groups_df, las = 1, main = "Species richness by group", ylab = "Number of species", xlab = "Groups") 
points(c(1:2), ntip_by_groups_colias, pch = 19, col = "red") 
legend("toprigh", bty = "n", legend = "Empirical values", pch = 19, col = "red", cex = 1)

median(ntip_by_groups_df$a)
median(ntip_by_groups_df$z)

library(coda)
divA <- as.mcmc(ntip_by_groups_df$a)
divZ <- as.mcmc(ntip_by_groups_df$z)

HPDinterval(divA)
HPDinterval(divZ)


# Tree imbalance

library(treebalance)

avgLeafDepI <- avgLeafDepI(tree) 
avgLeafDepI_posteriors <- sapply(posterior_trees, avgLeafDepI) 
hist(avgLeafDepI_posteriors, 
     main = "Distribution of average leaf depth index",
     xlab = "Average leaf depth", las = 1) 
abline(v = avgLeafDepI, col = "red")


# LTT 

library(phytools)

ltt_colias <- ltt(tree, plot = F) 
ltt_colias_df <- data.frame(times = round(ltt_colias$times, 4), ltt = ltt_colias$ltt)

posterior_trees_mp <- as.multiPhylo(posterior_trees[[1]]) 

for(i in 2:length(posterior_trees)){ 
  posterior_trees_mp[[i]] <- posterior_trees[[i]] 
} 

ltt95_CI <- ltt95(posterior_trees_mp, log = T, las = 1) 
lines(ltt_colias_df$times, ltt_colias_df$ltt, type = "s", col = "red") 
mtext(text = "Colias", side = 3, line = 1)

ltt95_CI_df <- as.data.frame(ltt95_CI[,c("time", "low(lineages)", "high(lineages)")])

ltt95_CI_df$time <- round(ltt95_CI_df$time, 4)

points_in_colias <- c() 
for(i in 1:nrow(ltt_colias_df)){
  int_max <- sort(ltt95_CI_df$time[ ltt95_CI_df$time >= ltt_colias_df$times[i]][1]) 
  int_min <- sort(ltt95_CI_df$time[ltt95_CI_df$time <= ltt_colias_df$times[i]],
                  decreasing = T)[1] 
  ltt_min <- ltt95_CI_df$`low(lineages)`[ltt95_CI_df$time == int_min] 
  ltt_max <- ltt95_CI_df$`high(lineages)`[ltt95_CI_df$time == int_max]
  points_in_colias[i] <- ifelse(ltt_min <= ltt_colias_df$ltt[i] & 
                                   ltt_colias_df$ltt[i] <= ltt_max, T, F)
} 

legend("topleft", legend = c("95% of the distribution of simulated trees around the median", "Median of simulated trees", "Empirical data"), lty = c(3,1,1), lwd = c(1,2,1), col = c("black","black","red"), bty = "n", cex = 0.8)

sum(points_in_colias)/length(points_in_colias)
# 98.5% of the values of the empirical LTT fall in the 95% confidence interval




# _ggplot ----
# library(tidyverse)
# 
# ntips <- bind_rows(ntip_by_groups_df, 
#                    tibble(a = ntip_by_groups_colias[1], z = ntip_by_groups_colias[2])) %>% 
#   mutate(data = c(rep("simulated", nrow(ntip_by_groups_df)), "empirical"),
#          id = 1:(nrow(ntip_by_groups_df)+1)) %>% 
#   rename(clade = a, backbone = z) %>% 
#   pivot_longer(1:2, names_to = "clade", values_to = "ntips")
# 
# ggplot(ntips) +
#   geom_violin(aes(clade, ntips, fill = data, col = data), alpha = 0.6) +
#   geom_point(aes(clade, ntips), data = filter(ntips, data == 'empirical')) +
#   theme_bw()
# 
# ggplot(ntips) +
#   geom_boxplot(aes(data, ntips), alpha = 0.6, data = filter(ntips, data == 'simulated')) +
#   geom_hline(aes(yintercept = ntips), col = "red", data = filter(ntips, data == 'empirical')) +
#   facet_grid(rows = vars(clade), scales = "free") +
#   theme_bw()




#--


#### Test model order ----------------------

shifts2 <- shift.estimates(phylo = tree,
                          data = groups,
                          sampling.fractions = sfrac,
                          comb.shift = node_comb,
                          # models = c("BCST","BVAR",
                          #            "BCST_DCST","BVAR_DVAR",
                          #            "BVAR_DCST","BCST_DVAR"),
                          models = c("BCST","BCST_DCST",
                                     "BVAR","BVAR_DCST",
                                     "BCST_DVAR","BVAR_DVAR"),
                          backbone.option = "crown.shift",
                          multi.backbone = "all",
                          Ncores = 4)

head(shifts2$total)
shifts2$backbones$`75/`

plot.phylo.comb(phylo = tree, 
                data = groups, 
                sampling.fractions = sfrac,
                shift.res = shifts2,
                combi = 1, 
                edge.width = 2, 
                label.offset = 0.3, 
                lad = F, 
                cex = 0.4)
axisPhylo()

rates2 <- div.rates(phylo = tree,
                   shift.res = shifts2,
                   combi = 1, part = "all")
rates2

paleodiversity2 <- paleodiv(phylo = tree,
                           data = groups,
                           sampling.fractions = sfrac,
                           shift.res = shifts2)
paleodiversity2

# compare

identical(sfrac, sfrac2)
identical(node_comb, node_comb2)
identical(shifts, shifts2)
identical(rates, rates2)



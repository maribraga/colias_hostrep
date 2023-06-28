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



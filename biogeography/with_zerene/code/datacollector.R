###########################
########## SETUP ##########
###########################

library(tools)

data_dir <- "../../output/sim/non_rj/sim_data"
rj_dir <- "../../output/sim/rj/rj_data"
geo_dir <- "../../data/sim/geo"


lio_dir <- "./biogeography/data_ignore/for_multifig/features"   #../../data/emp/features
out_dir <- "./biogeography/server_ignore/output_50k/output2"    #../../output/processed_data 

files <- list.files(path=data_dir, pattern=paste("*", ".summary.csv", sep=""), full.names=FALSE, recursive=FALSE)
rj_files <- list.files(path=rj_dir, pattern=paste("*", ".summary.csv", sep=""), full.names=FALSE, recursive=FALSE)

param_list <- c("r_b","r_d","r_e","r_w","sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]")
rj_param_list <- c("r_b","r_d","r_e","r_w","sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]","rj_sigma_w[1]","rj_sigma_w[2]","rj_sigma_e[1]","rj_sigma_e[2]","rj_sigma_b[1]","rj_sigma_b[2]","rj_sigma_d[1]","rj_sigma_d[2]","rj_phi_w[1]","rj_phi_w[2]","rj_phi_e[1]","rj_phi_e[2]","rj_phi_b[1]","rj_phi_b[2]","rj_phi_d[1]","rj_phi_d[2]")
geo_features <- c("qw_area","qw_altitude","qb_distance","qb_altitude")

######################################
########## SIM OUTPUT FILES ##########
######################################

# OUTPUT FORMAT:
# SIM,NREGIONS,GEONUM,TREESIZE,SUCCESS,PARAM1.MEAN,PARAM1.LOW80,PARAM1.HIGH80,PARAM1.TRUE,PARAM2.MEAN,PARAM2.LOW80,PARAM2.HIGH80,PARAM2.TRUE, ...

sim_list <- c()

for (n in c(3,5,7)) {
  for (g in 1:150) {
    sim <- paste(n,g,"medium",sep=".")
    sim_list <- c(sim_list, sim)
  }
}
for (g in 1:150) {
  sim <- paste(5,g,"xsmall",sep=".")
  sim_list <- c(sim_list, sim)
  sim <- paste(5,g,"small",sep=".")
  sim_list <- c(sim_list, sim)
  sim <- paste(5,g,"large",sep=".")
  sim_list <- c(sim_list, sim)
}

nrows <- length(sim_list)
blank_col <- rep(NA, nrows)
dataframe <- data.frame(sim_list,row.names=NULL)
rj_dataframe <- data.frame(sim_list,row.names=NULL)

nregions <- blank_col
geonum <- blank_col
treesize <- blank_col
success <- blank_col
dataframe <- cbind(dataframe,nregions,geonum,treesize,success)
rj_dataframe <- cbind(rj_dataframe,nregions,geonum,treesize,success)

for (param in param_list) {
  dataframe <- cbind(dataframe, blank_col)
  colnames(dataframe)[ncol(dataframe)] <- paste(param, ".mean", sep="")
  dataframe <- cbind(dataframe, blank_col)
  colnames(dataframe)[ncol(dataframe)] <- paste(param, ".low80", sep="")
  dataframe <- cbind(dataframe, blank_col)
  colnames(dataframe)[ncol(dataframe)] <- paste(param, ".high80", sep="")
  dataframe <- cbind(dataframe, blank_col)
  colnames(dataframe)[ncol(dataframe)] <- paste(param, ".true", sep="")
}

for (param in rj_param_list) {
  rj_dataframe <- cbind(rj_dataframe, blank_col)
  colnames(rj_dataframe)[ncol(rj_dataframe)] <- paste(param, ".mean", sep="")
  rj_dataframe <- cbind(rj_dataframe, blank_col)
  colnames(rj_dataframe)[ncol(rj_dataframe)] <- paste(param, ".low80", sep="")
  rj_dataframe <- cbind(rj_dataframe, blank_col)
  colnames(rj_dataframe)[ncol(rj_dataframe)] <- paste(param, ".high80", sep="")
  rj_dataframe <- cbind(rj_dataframe, blank_col)
  colnames(rj_dataframe)[ncol(rj_dataframe)] <- paste(param, ".true", sep="")
}

for (sim in sim_list) {
  index <- which(dataframe$sim_list==sim)
  filename <- paste(sim,".summary.csv",sep="")
  specs <- strsplit(sim,".",fixed=TRUE)
  n_regions <- specs[[1]][1]
  geo_number <- specs[[1]][2]
  tree_size <- specs[[1]][3]
  # NON-RJ
  dataframe$nregions[index] <- n_regions
  dataframe$geonum[index] <- geo_number
  dataframe$treesize[index] <- tree_size
  success <- filename %in% files
  dataframe$success[index] <- success
  if (success) {
    data <- read.csv(file=paste(data_dir, "/", sim, ".summary.csv", sep=""), sep = ",", header=TRUE)
    for (param in param_list) {
      param_ind <- which(data$parameter==param)
      mean_ind <- which(colnames(dataframe)==paste(param,".mean",sep=""))
      mean <- data[param_ind,2]
      dataframe[index,mean_ind] <- mean
      low_ind <- which(colnames(dataframe)==paste(param,".low80",sep=""))
      low <- data[param_ind,3]
      dataframe[index,low_ind] <- low
      high_ind <- which(colnames(dataframe)==paste(param,".high80",sep=""))
      high <- data[param_ind,4]
      dataframe[index,high_ind] <- high
      true_ind <- which(colnames(dataframe)==paste(param,".true",sep=""))
      true <- data[param_ind,5]
      dataframe[index,true_ind] <- true
    }
  }
  # RJ
  rj_dataframe$nregions[index] <- n_regions
  rj_dataframe$geonum[index] <- geo_number
  rj_dataframe$treesize[index] <- tree_size
  success <- filename %in% rj_files
  rj_dataframe$success[index] <- success
  if (success) {
    data <- read.csv(file=paste(rj_dir, "/", sim, ".summary.csv", sep=""), sep = ",", header=TRUE)
    for (param in rj_param_list) {
      param_ind <- which(data$parameter==param)
      mean_ind <- which(colnames(rj_dataframe)==paste(param,".mean",sep=""))
      mean <- data[param_ind,2]
      rj_dataframe[index,mean_ind] <- mean
      low_ind <- which(colnames(rj_dataframe)==paste(param,".low80",sep=""))
      low <- data[param_ind,3]
      rj_dataframe[index,low_ind] <- low
      high_ind <- which(colnames(rj_dataframe)==paste(param,".high80",sep=""))
      high <- data[param_ind,4]
      rj_dataframe[index,high_ind] <- high
      true_ind <- which(colnames(rj_dataframe)==paste(param,".true",sep=""))
      true <- data[param_ind,5]
      rj_dataframe[index,true_ind] <- true
    }
  }
}

##############################
########## GEO FILE ##########
##############################

geo_list <- c()

for (n in c(3,5,7)) {
  for (g in 1:150){
    geo <- paste(n,".",g,sep="")
    geo_list <- c(geo_list,geo)
  }
}

nrows <- length(geo_list)
blank_col <- rep(NA, nrows)
nregions <- blank_col
geonum <- blank_col
geo_dataframe <- data.frame(geo_list,nregions,geonum,row.names=NULL)

for (feature in geo_features) {
  geo_dataframe <- cbind(geo_dataframe, blank_col)
  colnames(geo_dataframe)[ncol(geo_dataframe)] <- paste(feature, ".variance", sep="")
}

for (geo in geo_list) {
  index <- which(geo_dataframe$geo_list==geo)
  specs <- strsplit(geo,".",fixed=TRUE)
  n_regions <- specs[[1]][1]
  geo_number <- specs[[1]][2]
  geo_dataframe$nregions[index] <- n_regions
  geo_dataframe$geonum[index] <- geo_number

  qw_area <- read.csv(paste(geo_dir,"/",n_regions,".",geo_number,".","qw_area.csv", sep=""),header = TRUE)
  qw_areas <- as.vector(t(qw_area))
  qw_area_gm <- exp(mean(log(abs(qw_areas))))
  qw_area_norm <- qw_areas/qw_area_gm
  qw_area_var <- var(qw_area_norm)
  geo_dataframe$qw_area.variance[index] <- qw_area_var

  qw_altitude <- read.csv(paste(geo_dir,"/",n_regions,".",geo_number,".","qw_altitude.csv", sep=""),header = TRUE)
  qw_altitudes <- as.vector(t(qw_altitude))
  qw_altitude_gm <- exp(mean(log(abs(qw_altitudes))))
  qw_altitude_norm <- qw_altitudes/qw_altitude_gm
  qw_altitude_var <- var(qw_altitude_norm)
  geo_dataframe$qw_altitude.variance[index] <- qw_altitude_var

  qb_distance <- read.csv(paste(geo_dir,"/",n_regions,".",geo_number,".","qb_distance.csv", sep=""),header = TRUE)
  qb_distances <- c()
  for (i in 1:n_regions) {
    for (j in 1:n_regions) {
      if (i != j) {
        qb_distances <- c(qb_distances, qb_distance[i,j])
      }
    }
  }
  qb_distance_gm <- exp(mean(log(abs(qb_distances))))
  qb_distance_norm <- abs(qb_distances)/qb_distance_gm
  qb_distance_var <- var(qb_distance_norm)
  geo_dataframe$qb_distance.variance[index] <- qb_distance_var

  qb_altitude <- read.csv(paste(geo_dir,"/",n_regions,".",geo_number,".","qb_altitude.csv", sep=""),header = TRUE)
  qb_altitudes <- c()
  for (i in 1:n_regions) {
    for (j in 1:n_regions) {
      if (i!=j) {
        qb_altitudes <- c(qb_altitudes, qb_altitude[i,j])
      }
    }
  }
  qb_altitude_gm <- exp(mean(log(abs(qb_altitudes))))
  qb_altitude_norm <- abs(qb_altitudes)/qb_altitude_gm
  qb_altitude_var <- var(qb_altitude_norm)
  geo_dataframe$qb_altitude.variance[index] <- qb_altitude_var
}

##############################
########## LIO FILE ##########
##############################

nregions <- c(7)

qw_area <- read.csv(paste(lio_dir,"/areas/areas.csv",sep=""),header=TRUE)
qw_areas <- as.vector(t(qw_area))
qw_area_gm <- exp(mean(log(abs(qw_areas))))
qw_area_norm <- qw_areas/qw_area_gm
qw_area_var <- var(qw_area_norm)
qw_area.variance <- qw_area_var

qw_altitude <- read.csv(paste(lio_dir,"/altitudes/mean.csv",sep=""),header=TRUE)
qw_altitudes <- as.vector(t(qw_altitude))
qw_altitude_gm <- exp(mean(log(abs(qw_altitudes))))
qw_altitude_norm <- qw_altitudes/qw_altitude_gm
qw_altitude_var <- var(qw_altitude_norm)
qw_altitude.variance <- qw_altitude_var

qb_distance <- read.csv(paste(lio_dir,"/distances/mean.csv",sep=""),header=TRUE)
qb_distances <- c()
for (i in 1:nregions) {
  for (j in 1:nregions) {
    if (i != j) {
      qb_distances <- c(qb_distances, qb_distance[i,j])
    }
  }
}
qb_distance_gm <- exp(mean(log(abs(qb_distances))))
qb_distance_norm <- abs(qb_distances)/qb_distance_gm
qb_distance_var <- var(qb_distance_norm)
qb_distance.variance <- qb_distance_var

qb_altitude <- read.csv(paste(lio_dir,"/altitudes/mean_diff.csv",sep=""),header=TRUE)
qb_altitudes <- c()
for (i in 1:nregions) {
  for (j in 1:nregions) {
    if (i!=j) {
      qb_altitudes <- c(qb_altitudes, qb_altitude[i,j])
    }
  }
}
qb_altitude_gm <- exp(mean(log(abs(qb_altitudes))))
qb_altitude_norm <- abs(qb_altitudes)/qb_altitude_gm
qb_altitude_var <- var(qb_altitude_norm)
qb_altitude.variance <- qb_altitude_var

lio_dataframe <- data.frame(nregions,qw_area.variance,qw_altitude.variance,qb_distance.variance,qb_altitude.variance,row.names=NULL)

##################################
########## SAVING FILES ##########
##################################

write.csv(dataframe,file=paste(out_dir, "/", "regular_sim.csv", sep=""),row.names=FALSE,quote=FALSE)
write.csv(rj_dataframe,file=paste(out_dir, "/", "rj_sim.csv", sep=""),row.names=FALSE,quote=FALSE)
write.csv(geo_dataframe,file=paste(out_dir, "/", "geo_variances.csv", sep=""),row.names=FALSE,quote=FALSE)
write.csv(lio_dataframe,file=paste(out_dir, "/", "lio_variances.csv", sep=""),row.names=FALSE,quote=FALSE)

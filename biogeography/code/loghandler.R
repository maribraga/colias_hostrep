library(coda)

process_logfile <- function(data, param) {
  param2 <- gsub("\\[|\\]", ".", param)
  param_data <- data[param2]
  param_data <- param_data[,1]
  # Removing burnin
  param_data <- param_data[101:1001]
  # Removing zeroes
  param_data <- param_data[!param_data == 0]
  # Getting mean
  mean <- mean(param_data)
  # Getting HPD
  coda <- HPDinterval(as.mcmc(param_data),prob=.8)
  low <- coda[1,1]
  high <- coda[1,2]
  return(c(param,mean,low,high))
}

log_dir <- "../../output/sim/non_rj/sim_output"
true_dir <- "../../data/sim/history"
summary_dir <- "../../output/sim/non_rj/sim_data"

files <- list.files(path=log_dir, pattern=paste("*", ".N.model.log", sep=""), full.names=TRUE, recursive=FALSE)
params <- c("r_b","r_d","r_e","r_w","sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]")

num <- 0
for (file in files) {
  num <- num + 1
  cat("\r",round(num/length(files)*100,0),"%")
  results <- data.frame(matrix(ncol=5,nrow=0))
  specs <- strsplit(basename(file),".",fixed=TRUE)
  n_regions <- specs[[1]][1]
  geo_number <- specs[[1]][2]
  tree_size <- specs[[1]][3]
  sim <- paste(n_regions, geo_number, tree_size, sep=".")
  data <- read.csv(file=paste(log_dir, "/", sim, ".N.model.log", sep=""), sep = "\t", header=TRUE)
  length <- length(data[,1])
  if (length == 1000 + 1) {
    for (param in params) {
      statistics <- process_logfile(data, param)
      true_data <- read.csv(file=paste(true_dir, "/", sim, ".param.txt", sep=""), sep = ",", header=TRUE)
      ind <- which(true_data$parameter == param)
      val <- true_data[ind,2]
      statistics <- c(statistics, val)
      results <- rbind(results, statistics)
    }
    colnames(results) <- c("parameter", "mean", "80hpd_low", "80hpd_high", "true")
    write.csv(results,file=paste(summary_dir, "/", sim, ".summary.csv", sep=""),row.names=FALSE,quote=FALSE)
  }
}

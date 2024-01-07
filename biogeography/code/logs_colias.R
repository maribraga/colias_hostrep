library(coda)
library(tidyverse)

process_logfile <- function(data, param) {
  param2 <- gsub("\\[|\\]", ".", param)
  param_data <- data[param2]
  param_data <- param_data[,1]
  # Removing burnin
  burnin <- floor(length(param_data)/10)
  param_data <- param_data[(burnin+1):length(param_data)]
  # Removing zeroes
  param_data <- param_data[!param_data == 0]
  # Getting mean
  mean <- mean(param_data)
  # Getting HPD
  coda <- HPDinterval(as.mcmc(param_data),prob=.95)
  low <- coda[1,1]
  high <- coda[1,2]
  
  return(data.frame(param, mean, low, high))
}

#params <- c("r_b","r_d","r_e","r_w","sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]")
params <- c("phi_b[1]","phi_b[2]","phi_b[3]","phi_d[1]","phi_d[2]","phi_d[3]","phi_e[1]","phi_e[2]","phi_e[3]","phi_w[1]","phi_w[2]","phi_w[3]","r_e[1]","r_e[2]","r_e[3]","r_e[4]","r_e[5]","r_e[6]","r_e[7]","rho_b","rho_d","rho_e","rho_w","sigma_b[1]","sigma_b[2]","sigma_b[3]","sigma_d[1]","sigma_d[2]","sigma_d[3]","sigma_e[1]","sigma_e[2]","sigma_e[3]","sigma_w[1]","sigma_w[2]","sigma_w[3]")

data <- read.csv(file="/Users/mari/repos/colias_hostrep/biogeography/server_ignore/output_50k/output2/model.log", sep = "\t", header=TRUE)
length <- length(data[,1])

results <- data.frame(matrix(ncol=4,nrow=0))

for (param in params) {
  statistics <- process_logfile(data, param)
  results <- rbind(results, statistics)
}
colnames(results) <- c("parameter", "mean", "hpd95_low", "hpd95_high")
write.csv(results,file="biogeography/result_summary.csv",row.names=FALSE,quote=FALSE)


ggplot(results) +
  geom_point(aes(parameter, mean)) +
  geom_linerange(aes(x = parameter, ymin = hpd95_low, ymax = hpd95_high)) +
  geom_hline(yintercept = 0.0, col = "red") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 90))
  

params_subset <- c("r_e[1]","r_e[2]","r_e[3]","r_e[4]","r_e[5]","r_e[6]","r_e[7]","rho_b","rho_d","rho_e","rho_w")

results %>% 
  filter(parameter %in% params_subset)

results %>% 
  filter(str_detect(parameter, "r_e")) %>% 
  ggplot() +
  geom_point(aes(parameter, mean)) +
  geom_linerange(aes(x = parameter, ymin = hpd95_low, ymax = hpd95_high)) +
  #geom_hline(yintercept = 0.0, col = "red") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 90))

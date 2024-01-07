###########################
########## SETUP ##########
###########################

# Packages
library(tools)
library(ggplot2)
library(patchwork)
library(stringr)
library(gt)
library(webshot2)
library(coda)
library(ggridges)
library(RevGadgets)
library(ggtree)
library(maptools)
library(rgeos)
library(dplyr)
library(igraph)
library(sf)
library(rnaturalearth)
library(RColorBrewer)
library(RevGadgets)
library(Polychrome)

########################################
########## PLOTTING FUNCTIONS ##########
########################################

make_varplot <- function () {
  breaks <- c(1e-02,1e-01,1e+00,1e+01,1e+02)
  labels <- c(".01",".1","1","10","100")
  limits=c(1e-2,1e+2)

  qw_area.qw_altitude <-
    ggplot(data=geo_variances_data,aes(x=qw_area.variance,y=qw_altitude.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qw_area.variance[1],y=lio_variances_data$qw_altitude.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Area}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  qw_area.qb_distance <-
    ggplot(data=geo_variances_data,aes(x=qw_area.variance,y=qb_distance.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qw_area.variance[1],y=lio_variances_data$qb_distance.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Distance}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  qw_area.qb_altitude <-
    ggplot(data=geo_variances_data,aes(x=qw_area.variance,y=qb_altitude.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qw_area.variance[1],y=lio_variances_data$qb_altitude.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Area}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  qw_altitude.qb_distance <-
    ggplot(data=geo_variances_data,aes(x=qw_altitude.variance,y=qb_distance.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qw_altitude.variance[1],y=lio_variances_data$qb_distance.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Distance}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  qw_altitude.qb_altitude <-
    ggplot(data=geo_variances_data,aes(x=qw_altitude.variance,y=qb_altitude.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qw_altitude.variance[1],y=lio_variances_data$qb_altitude.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(w)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  qb_distance.qb_altitude <-
    ggplot(data=geo_variances_data,aes(x=qb_distance.variance,y=qb_altitude.variance)) +
    geom_point(size=.5) +
    geom_point(aes(x=lio_variances_data$qb_distance.variance[1],y=lio_variances_data$qb_altitude.variance[1]),color="red",size=2) +
    scale_x_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Distance}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    scale_y_continuous(name=bquote("Variance of" ~ bold(b)[q]^{Altitude}),breaks=breaks,labels=labels,limits=limits,trans="log10") +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_text(angle=90,vjust=0,hjust=.5))

  variance_plot <-
    qw_area.qw_altitude + qw_area.qb_distance + qw_area.qb_altitude + qw_altitude.qb_distance + qw_altitude.qb_altitude + qb_distance.qb_altitude +
    plot_annotation(tag_levels="a") &
    theme(plot.tag=element_text(face="bold"))

  return(variance_plot)
}

make_covplots <- function () {
  coverages <- data.frame(matrix(ncol=7,nrow=0))
  m5.all <- data.frame(matrix(ncol=6,nrow=0))
  for (param in params) {
    param2 <- gsub("\\[|\\]", ".", param)
    low_colname <- paste(param2,".low80",sep="")
    true_colname <- paste(param2,".true",sep="")
    high_colname <- paste(param2,".high80",sep="")
    est_colname <- paste(param2,".mean",sep="")
    # M3
    m3.low <- m3[,low_colname]
    m3.true <- m3[,true_colname]
    m3.high <- m3[,high_colname]
    m3.est <- m3[,est_colname]
    m3.suc <- rep(FALSE,nrow(m3))
    m3.cov <- data.frame(m3.low,m3.true,m3.high,m3.est,m3.suc)
    for (i in 1:nrow(m3.cov)) {
      if (m3.cov[i,1] < m3.cov[i,2] & m3.cov[i,2] < m3.cov[i,3]) {
        m3.cov[i,5] <- TRUE
      }
    }
    m3.c <- length(which(m3.cov$m3.suc==TRUE))/nrow(m3.cov)
    # x5
    x5.low <- x5[,low_colname]
    x5.true <- x5[,true_colname]
    x5.high <- x5[,high_colname]
    x5.est <- x5[,est_colname]
    x5.suc <- rep(FALSE,nrow(x5))
    x5.cov <- data.frame(x5.low,x5.true,x5.high,x5.est,x5.suc)
    for (i in 1:nrow(x5.cov)) {
      if (x5.cov[i,1] < x5.cov[i,2] & x5.cov[i,2] < x5.cov[i,3]) {
        x5.cov[i,5] <- TRUE
      }
    }
    x5.c <- length(which(x5.cov$x5.suc==TRUE))/nrow(x5.cov)
    # s5
    s5.low <- s5[,low_colname]
    s5.true <- s5[,true_colname]
    s5.high <- s5[,high_colname]
    s5.est <- s5[,est_colname]
    s5.suc <- rep(FALSE,nrow(s5))
    s5.cov <- data.frame(s5.low,s5.true,s5.high,s5.est,s5.suc)
    for (i in 1:nrow(s5.cov)) {
      if (s5.cov[i,1] < s5.cov[i,2] & s5.cov[i,2] < s5.cov[i,3]) {
        s5.cov[i,5] <- TRUE
      }
    }
    s5.c <- length(which(s5.cov$s5.suc==TRUE))/nrow(s5.cov)
    # m5
    m5.low <- m5[,low_colname]
    m5.true <- m5[,true_colname]
    m5.high <- m5[,high_colname]
    m5.est <- m5[,est_colname]
    m5.suc <- rep(FALSE,nrow(m5))
    m5.cov <- data.frame(m5.low,m5.true,m5.high,m5.est,m5.suc)
    for (i in 1:nrow(m5.cov)) {
      if (m5.cov[i,1] < m5.cov[i,2] & m5.cov[i,2] < m5.cov[i,3]) {
        m5.cov[i,5] <- TRUE
      }
    }
    m5.c <- length(which(m5.cov$m5.suc==TRUE))/nrow(m5.cov)
    m5.param <- rep(param,nrow(m5))
    m5.param.all <- cbind(m5.cov,m5.param)
    m5.all <- rbind(m5.all,m5.param.all)
    # l5
    l5.low <- l5[,low_colname]
    l5.true <- l5[,true_colname]
    l5.high <- l5[,high_colname]
    l5.est <- l5[,est_colname]
    l5.suc <- rep(FALSE,nrow(l5))
    l5.cov <- data.frame(l5.low,l5.true,l5.high,l5.est,l5.suc)
    for (i in 1:nrow(l5.cov)) {
      if (l5.cov[i,1] < l5.cov[i,2] & l5.cov[i,2] < l5.cov[i,3]) {
        l5.cov[i,5] <- TRUE
      }
    }
    l5.c <- length(which(l5.cov$l5.suc==TRUE))/nrow(l5.cov)
    # m7
    m7.low <- m7[,low_colname]
    m7.true <- m7[,true_colname]
    m7.high <- m7[,high_colname]
    m7.est <- m7[,est_colname]
    m7.suc <- rep(FALSE,nrow(m7))
    m7.cov <- data.frame(m7.low,m7.true,m7.high,m7.est,m7.suc)
    for (i in 1:nrow(m7.cov)) {
      if (m7.cov[i,1] < m7.cov[i,2] & m7.cov[i,2] < m7.cov[i,3]) {
        m7.cov[i,5] <- TRUE
      }
    }
    m7.c <- length(which(m7.cov$m7.suc==TRUE))/nrow(m7.cov)
    coverage <- c(param_names[which(params==param)], round(m3.c*100,2), round(x5.c*100,2), round(s5.c*100,2), round(m5.c*100,2), round(l5.c*100,2), round(m7.c*100,2))
    coverages <- rbind(coverages, coverage)
  }
  colnames(coverages) <- c("parameter","m3","x5","s5","m5","l5","m7")
  colnames(m5.all) <- c("low","true","high","est","suc","param")

  plot_param_rhol <- function(data, name) {
    min <- 0
    max <- 6
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,8,1,"pt"))
    return(plot)
  }

  plot_param_rho <- function(data, name) {
    min <- 0
    max <- 6
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,8,1,"pt"),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    return(plot)
  }

  plot_param_l <- function(data, name) {
    min <- -3
    max <- 3
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,1,1,"pt"),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    return(plot)
  }

  plot_param_b <- function(data, name) {
    min <- -3
    max <- 3
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,1,1,"pt"),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    return(plot)
  }

  plot_param_lb <- function(data, name) {
    min <- -3
    max <- 3
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,1,1,"pt"))
    return(plot)
  }

  plot_param <- function(data, name) {
    min <- -3
    max <- 3
    cov <- round(length(which(data$suc==TRUE))/length(data$suc)*100)
    plot <- ggplot(data, aes(x=true, y=est)) +
      geom_errorbar(aes(ymin=ifelse(low>min,low,min),ymax=ifelse(high<max,high,max),col=suc),size=.25) +
      geom_point(size=.5) +
      geom_abline(slope=1,size=.25) +
      labs(x=NULL,y=NULL) +
      annotate(geom="text",label=name,x=min,y=max,hjust="inward",vjust="inward") +
      annotate(geom="text",label=paste(cov,"%",sep=""),x=max,y=min,hjust="inward",vjust="inward") +
      lims(x=c(min,max),y=c(min,max)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(1,1,1,1,"pt"),axis.text=element_blank(),axis.ticks=element_blank())
    return(plot)
  }

  rw <- m5.all[which(m5.all$param=="r_w"),]
  rw <- plot_param_rhol(rw, bquote(rho[w]))
  re <- m5.all[which(m5.all$param=="r_e"),]
  re <- plot_param_rho(re, bquote(rho[e]))
  rb <- m5.all[which(m5.all$param=="r_b"),]
  rb <- plot_param_rho(rb, bquote(rho[b]))
  rd <- m5.all[which(m5.all$param=="r_d"),]
  rd <- plot_param_rho(rd, bquote(rho[d]))

  sigmaw1 <- m5.all[which(m5.all$param=="sigma_w[1]"),]
  sigmaw1 <- plot_param_l(sigmaw1, bquote(sigma[w]^{Area}))
  sigmaw2 <- m5.all[which(m5.all$param=="sigma_w[2]"),]
  sigmaw2 <- plot_param(sigmaw2, bquote(sigma[w]^{Altitude}))
  sigmae1 <- m5.all[which(m5.all$param=="sigma_e[1]"),]
  sigmae1 <- plot_param(sigmae1, bquote(sigma[e]^{Area}))
  sigmae2 <- m5.all[which(m5.all$param=="sigma_e[2]"),]
  sigmae2 <- plot_param(sigmae2, bquote(sigma[e]^{Altitude}))
  sigmab1 <- m5.all[which(m5.all$param=="sigma_b[1]"),]
  sigmab1 <- plot_param_l(sigmab1, bquote(sigma[b]^{Distance}))
  sigmab2 <- m5.all[which(m5.all$param=="sigma_b[2]"),]
  sigmab2 <- plot_param(sigmab2, bquote(sigma[b]^{Altitude}))
  sigmad1 <- m5.all[which(m5.all$param=="sigma_d[1]"),]
  sigmad1 <- plot_param(sigmad1, bquote(sigma[d]^{Distance}))
  sigmad2 <- m5.all[which(m5.all$param=="sigma_d[2]"),]
  sigmad2 <- plot_param(sigmad2, bquote(sigma[d]^{Altitude}))
  phiw1 <- m5.all[which(m5.all$param=="phi_w[1]"),]
  phiw1 <- plot_param_l(phiw1, bquote(phi[w]^{Area}))
  phiw2 <- m5.all[which(m5.all$param=="phi_w[2]"),]
  phiw2 <- plot_param(phiw2, bquote(phi[w]^{Altitude}))
  phie1 <- m5.all[which(m5.all$param=="phi_e[1]"),]
  phie1 <- plot_param(phie1, bquote(phi[e]^{Area}))
  phie2 <- m5.all[which(m5.all$param=="phi_e[2]"),]
  phie2 <- plot_param(phie2, bquote(phi[e]^{Altitude}))
  phib1 <- m5.all[which(m5.all$param=="phi_b[1]"),]
  phib1 <- plot_param_lb(phib1, bquote(phi[b]^{Distance}))
  phib2 <- m5.all[which(m5.all$param=="phi_b[2]"),]
  phib2 <- plot_param_b(phib2, bquote(phi[b]^{Altitude}))
  phid1 <- m5.all[which(m5.all$param=="phi_d[1]"),]
  phid1 <- plot_param_b(phid1, bquote(phi[d]^{Distance}))
  phid2 <- m5.all[which(m5.all$param=="phi_d[2]"),]
  phid2 <- plot_param_b(phid2, bquote(phi[d]^{Altitude}))

  axis_plot <- ggplot() +
    labs(x="True (Simulated) Parameter Value",y="Posterior Mean and 80% HPD Interval") +
    theme_classic() +
    theme(aspect.ratio=5/4,axis.title=element_text(hjust=.5,size=16),line=element_blank())
  inner_plot <- rw + re + rb + rd + sigmaw1 + sigmaw2 + sigmae1 + sigmae2 + sigmab1 + sigmab2 + sigmad1 + sigmad2 + phiw1 + phiw2 + phie1 + phie2 + phib1 + phib2 + phid1 + phid2 +
    plot_layout(ncol=4)
  inner_plot
  coverage_plot <-
    axis_plot + inset_element(inner_plot,left=0,bottom=0,right=1,top=1)

  coverages_table <- gt(coverages) %>%
    cols_label(parameter="Parameter",m3="3M",x5="5X",s5="5S",m5="5M",l5="5L",m7="7M") %>%
    cols_align(align="center",columns=everything())

  return(list(plot=coverage_plot,table=coverages_table))
}

make_rjplot <- function () {
  m5.rj <- data.frame(matrix(ncol=6,nrow=0))
  for (param in rjparams) {
    param2 <- gsub("\\[|\\]", ".", param)
    param3 <- paste("rj_",param2,sep="")
    true_colname <- paste(param2,".true",sep="")
    est_colname <- paste(param3,".mean",sep="")
    # m5
    m5rj.true <- abs(m5rj[,true_colname])
    m5rj.est <- m5rj[,est_colname]
    m5.dat <- data.frame(m5rj.true,m5rj.est)
    m5.param <- rep(param,nrow(m5rj))
    m5.param.all <- cbind(m5.dat,m5.param)
    m5.rj <- rbind(m5.rj,m5.param.all)
  }
  colnames(m5.rj) <- c("true","rjprob","param")

  plot_param <- function(data, name) {
    min <- min(min(data$true),min(data$rjprob))
    max <- max(max(data$true),max(data$rjprob))
    plot <- ggplot(data, aes(x=true, y=rjprob)) +
      geom_point(size=.5) +
      geom_smooth(method='lm',size=.25) +
      annotate(geom="text",label=name,x=3,y=0,hjust="inward",vjust="inward") +
      labs(x=NULL,y=NULL) +
      lims(x=c(0,3),y=c(0,1)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(2,2,2,2,"pt"),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())
    return(plot)
  }

  plot_param_left <- function(data, name) {
    min <- min(min(data$true),min(data$rjprob))
    max <- max(max(data$true),max(data$rjprob))
    plot <- ggplot(data, aes(x=true, y=rjprob)) +
      geom_point(size=.5) +
      geom_smooth(method='lm',size=.25) +
      annotate(geom="text",label=name,x=3,y=0,hjust="inward",vjust="inward") +
      labs(x=NULL,y=NULL) +
      lims(x=c(0,3),y=c(0,1)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(2,2,2,2,"pt"),axis.ticks.x=element_blank(),axis.text.x=element_blank())
    return(plot)
  }

  plot_param_bottom <- function(data, name) {
    min <- min(min(data$true),min(data$rjprob))
    max <- max(max(data$true),max(data$rjprob))
    plot <- ggplot(data, aes(x=true, y=rjprob)) +
      geom_point(size=.5) +
      geom_smooth(method='lm',size=.25) +
      annotate(geom="text",label=name,x=3,y=0,hjust="inward",vjust="inward") +
      labs(x=NULL,y=NULL) +
      lims(x=c(0,3),y=c(0,1)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(2,2,2,2,"pt"),axis.ticks.y=element_blank(),axis.text.y=element_blank())
    return(plot)
  }

  plot_param_bottomleft <- function(data, name) {
    min <- min(min(data$true),min(data$rjprob))
    max <- max(max(data$true),max(data$rjprob))
    plot <- ggplot(data, aes(x=true, y=rjprob)) +
      geom_point(size=.5) +
      geom_smooth(method='lm',size=.25) +
      annotate(geom="text",label=name,x=3,y=0,hjust="inward",vjust="inward") +
      labs(x=NULL,y=NULL) +
      lims(x=c(0,3),y=c(0,1)) +
      theme_bw() +
      theme(aspect.ratio=1,panel.grid.minor=element_blank(),panel.grid.major=element_blank(),legend.position="none",plot.margin=margin(2,2,2,2,"pt"))
    return(plot)
  }

  sigmaw1 <- m5.rj[which(m5.rj$param=="sigma_w[1]"),]
  sigmaw1 <- plot_param_left(sigmaw1, bquote(sigma[w]^{Area}))
  sigmaw2 <- m5.rj[which(m5.rj$param=="sigma_w[2]"),]
  sigmaw2 <- plot_param(sigmaw2, bquote(sigma[w]^{Altitude}))
  sigmae1 <- m5.rj[which(m5.rj$param=="sigma_e[1]"),]
  sigmae1 <- plot_param(sigmae1, bquote(sigma[e]^{Area}))
  sigmae2 <- m5.rj[which(m5.rj$param=="sigma_e[2]"),]
  sigmae2 <- plot_param(sigmae2, bquote(sigma[e]^{Altitude}))
  sigmab1 <- m5.rj[which(m5.rj$param=="sigma_b[1]"),]
  sigmab1 <- plot_param_left(sigmab1, bquote(sigma[b]^{Distance}))
  sigmab2 <- m5.rj[which(m5.rj$param=="sigma_b[2]"),]
  sigmab2 <- plot_param(sigmab2, bquote(sigma[b]^{Altitude}))
  sigmad1 <- m5.rj[which(m5.rj$param=="sigma_d[1]"),]
  sigmad1 <- plot_param(sigmad1, bquote(sigma[d]^{Distance}))
  sigmad2 <- m5.rj[which(m5.rj$param=="sigma_d[2]"),]
  sigmad2 <- plot_param(sigmad2, bquote(sigma[d]^{Altitude}))
  phiw1 <- m5.rj[which(m5.rj$param=="phi_w[1]"),]
  phiw1 <- plot_param_left(phiw1, bquote(phi[w]^{Area}))
  phiw2 <- m5.rj[which(m5.rj$param=="phi_w[2]"),]
  phiw2 <- plot_param(phiw2, bquote(phi[w]^{Altitude}))
  phie1 <- m5.rj[which(m5.rj$param=="phi_e[1]"),]
  phie1 <- plot_param(phie1, bquote(phi[e]^{Area}))
  phie2 <- m5.rj[which(m5.rj$param=="phi_e[2]"),]
  phie2 <- plot_param(phie2, bquote(phi[e]^{Altitude}))
  phib1 <- m5.rj[which(m5.rj$param=="phi_b[1]"),]
  phib1 <- plot_param_bottomleft(phib1, bquote(phi[b]^{Distance}))
  phib2 <- m5.rj[which(m5.rj$param=="phi_b[2]"),]
  phib2 <- plot_param_bottom(phib2, bquote(phi[b]^{Altitude}))
  phid1 <- m5.rj[which(m5.rj$param=="phi_d[1]"),]
  phid1 <- plot_param_bottom(phid1, bquote(phi[d]^{Distance}))
  phid2 <- m5.rj[which(m5.rj$param=="phi_d[2]"),]
  phid2 <- plot_param_bottom(phid2, bquote(phi[d]^{Altitude}))

  axis_plot <- ggplot() +
    labs(x="True (Simulated) Effect Strength",y="RJMCMC Probability") +
    theme_classic() +
    theme(line=element_blank())
  inner_plot <- sigmaw1 + sigmaw2 + sigmae1 + sigmae2 + sigmab1 + sigmab2 + sigmad1 + sigmad2 + phiw1 + phiw2 + phie1 + phie2 + phib1 + phib2 + phid1 + phid2 +
    plot_layout(ncol=4)
  rjprob_plot <-
    axis_plot + inset_element(inner_plot,left=0,bottom=0,right=1,top=1)

  return(rjprob_plot)
}

make_postplots <- function () {
  mcmc_nonrj <- as.mcmc(data_nonrj)
  hpd_nonrj <- HPDinterval(mcmc_nonrj,prob=.95)
  mcmc_rj <- as.mcmc(data_rj)
  hpd_rj <- HPDinterval(mcmc_rj,prob=.95)

  nonrj_post_rw <- ggplot(data_nonrj, aes(x=r_w,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["r_w",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["r_w",2],linetype="dashed",linewidth=.25) +
    geom_density() +
    labs(y="Posterior Density", x=bquote("Value of " ~ rho[w])) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(8,8,8,8))

  nonrj_post_re <- ggplot(data_nonrj, aes(x=r_e,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["r_e",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["r_e",2],linetype="dashed",linewidth=.25) +
    geom_density() +
    labs(y=NULL, x=bquote("Value of " ~ rho[e])) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(8,8,8,8))

  nonrj_post_rb <- ggplot(data_nonrj, aes(x=r_b,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["r_b",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["r_b",2],linetype="dashed",linewidth=.25) +
    geom_density() +
    labs(y=NULL, x=bquote("Value of " ~ rho[b])) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(8,8,8,8))

  nonrj_post_rd <- ggplot(data_nonrj, aes(x=r_d,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["r_d",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["r_d",2],linetype="dashed",linewidth=.25) +
    geom_density() +
    labs(y=NULL, x=bquote("Value of " ~ rho[d])) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(8,8,8,8))

  baserates_plot <- nonrj_post_rw + nonrj_post_re + nonrj_post_rb + nonrj_post_rd +
    plot_layout(ncol=4)

  nonrj_post_sigmaw1 <- ggplot(data_nonrj, aes(x=sigma_w.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_w.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_w.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Area " ~ sigma[w]^{Area})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmaw1 <- ggplot(data_rj, aes(x=rj_sigma_w.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phiw1 <- ggplot(data_nonrj, aes(x=phi_w.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_w.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_w.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Area " ~ phi[w]^{Area})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phiw1 <- ggplot(data_rj, aes(x=rj_phi_w.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_sigmaw2 <- ggplot(data_nonrj, aes(x=sigma_w.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_w.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_w.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Altitude " ~ sigma[w]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmaw2 <- ggplot(data_rj, aes(x=rj_sigma_w.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phiw2 <- ggplot(data_nonrj, aes(x=phi_w.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_w.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_w.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Altitude " ~ phi[w]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phiw2 <- ggplot(data_rj, aes(x=rj_phi_w.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  w <- nonrj_post_sigmaw1 + rj_prob_sigmaw1 + nonrj_post_phiw1 + rj_prob_phiw1 + nonrj_post_sigmaw2 + rj_prob_sigmaw2 + nonrj_post_phiw2 + rj_prob_phiw2 +
    plot_layout(ncol=4) +
    plot_annotation(theme=theme(plot.title=element_text(size=24,hjust=.5))) &
    theme(text=element_text(size=8))

  nonrj_post_sigmae1 <- ggplot(data_nonrj, aes(x=sigma_e.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_e.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_e.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Area " ~ sigma[e]^{Area})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmae1 <- ggplot(data_rj, aes(x=rj_sigma_e.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phie1 <- ggplot(data_nonrj, aes(x=phi_e.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_e.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_e.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Area " ~ phi[e]^{Area})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phie1 <- ggplot(data_rj, aes(x=rj_phi_e.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_sigmae2 <- ggplot(data_nonrj, aes(x=sigma_e.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_e.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_e.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Altitude " ~ sigma[e]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmae2 <- ggplot(data_rj, aes(x=rj_sigma_e.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phie2 <- ggplot(data_nonrj, aes(x=phi_e.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_e.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_e.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Altitude " ~ phi[e]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phie2 <- ggplot(data_rj, aes(x=rj_phi_e.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  e <- nonrj_post_sigmae1 + rj_prob_sigmae1 + nonrj_post_phie1 + rj_prob_phie1 + nonrj_post_sigmae2 + rj_prob_sigmae2 + nonrj_post_phie2 + rj_prob_phie2 +
    plot_layout(ncol=4) +
    plot_annotation(theme=theme(plot.title=element_text(size=24,hjust=.5))) &
    theme(text=element_text(size=8))

  nonrj_post_sigmab1 <- ggplot(data_nonrj, aes(x=sigma_b.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_b.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_b.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Distance " ~ sigma[b]^{Distance})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmab1 <- ggplot(data_rj, aes(x=rj_sigma_b.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phib1 <- ggplot(data_nonrj, aes(x=phi_b.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_b.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_b.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Distance " ~ phi[b]^{Distance})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phib1 <- ggplot(data_rj, aes(x=rj_phi_b.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_sigmab2 <- ggplot(data_nonrj, aes(x=sigma_b.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_b.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_b.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Altitude " ~ sigma[b]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmab2 <- ggplot(data_rj, aes(x=rj_sigma_b.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phib2 <- ggplot(data_nonrj, aes(x=phi_b.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_b.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_b.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Altitude " ~ phi[b]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phib2 <- ggplot(data_rj, aes(x=rj_phi_b.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  b <- nonrj_post_sigmab1 + rj_prob_sigmab1 + nonrj_post_phib1 + rj_prob_phib1 + nonrj_post_sigmab2 + rj_prob_sigmab2 + nonrj_post_phib2 + rj_prob_phib2 +
    plot_layout(ncol=4) +
    plot_annotation(theme=theme(plot.title=element_text(size=24,hjust=.5))) &
    theme(text=element_text(size=8))

  nonrj_post_sigmad1 <- ggplot(data_nonrj, aes(x=sigma_d.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_d.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_d.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Distance " ~ sigma[d]^{Distance})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmad1 <- ggplot(data_rj, aes(x=rj_sigma_d.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phid1 <- ggplot(data_nonrj, aes(x=phi_d.1.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_d.1.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_d.1.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Distance " ~ phi[d]^{Distance})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phid1 <- ggplot(data_rj, aes(x=rj_phi_d.1.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_sigmad2 <- ggplot(data_nonrj, aes(x=sigma_d.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["sigma_d.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["sigma_d.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Categorical Altitude " ~ sigma[d]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_sigmad2 <- ggplot(data_rj, aes(x=rj_sigma_d.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  nonrj_post_phid2 <- ggplot(data_nonrj, aes(x=phi_d.2.,y=..scaled..)) +
    geom_vline(xintercept=hpd_nonrj["phi_d.2.",1],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=hpd_nonrj["phi_d.2.",2],linetype="dashed",linewidth=.25) +
    geom_vline(xintercept=0,linewidth=.5,color="red") +
    geom_density() +
    labs(x="Parameter Estimate",y="Posterior Density",title=bquote("Quantitative Altitude " ~ phi[d]^{Altitude})) +
    theme_bw() +
    theme(aspect.ratio=1,panel.grid.minor=element_blank(),plot.margin=margin(4,2,4,4))
  rj_prob_phid2 <- ggplot(data_rj, aes(x=rj_phi_d.2.)) +
    geom_bar(fill="gray") +
    lims(y=c(0,2001)) +
    scale_x_continuous(breaks=c(0,1)) +
    labs(x="RJMCMC",y=NULL) +
    annotate(geom="text",size=2.25,label="Effect",x=1,y=2001) +
    annotate(geom="text",size=2.25,label="Non-Effect",x=0,y=2001) +
    theme_bw() +
    theme(aspect.ratio=2,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.margin=margin(4,4,4,2))

  d <- nonrj_post_sigmad1 + rj_prob_sigmad1 + nonrj_post_phid1 + rj_prob_phid1 + nonrj_post_sigmad2 + rj_prob_sigmad2 + nonrj_post_phid2 + rj_prob_phid2 +
    plot_layout(ncol=4) +
    plot_annotation(theme=theme(plot.title=element_text(size=24,hjust=.5))) &
    theme(text=element_text(size=8))

  return(list(r=baserates_plot,w=w,e=e,b=b,d=d))
}

make_jointplot <- function () {
  joint <- function(param1, param2) {
    param1b <- paste("rj_",gsub("\\[|\\]", ".", param1),sep="")
    param2b <- paste("rj_",gsub("\\[|\\]", ".", param2),sep="")
    rj1 <- data_rj[,param1b]
    rj2 <- data_rj[,param2b]
    joints <- c(0,0,0,0)
    for (i in 1:length(rj1)) {
      if (rj1[i] == 1) {
        if (rj2[i] == 1) {joints[1] <- joints[1] + 1}
        else {joints[2] <- joints[2] + 1}
      }
      else {
        if (rj2[i] == 1) {joints[3] <- joints[3] + 1}
        else {joints[4] <- joints[4] + 1}
      }
    }
    joints <- round(joints/sum(joints),2)
    return(joints)
  }

  jointplot <- function(param1,param2,name) {
    joints <- joint(param1,param2)
    dataframe <- data.frame(matrix(c(0,1,2,0,1,2),ncol=2))
    plot <- ggplot(dataframe,aes(x=X1,y=X2)) +
      geom_rect(xmin=1,xmax=2,ymin=1,ymax=2,color="black",fill="blue",alpha=joints[1]) +
      geom_rect(xmin=1,xmax=2,ymin=0,ymax=1,color="black",fill="blue",alpha=joints[2]) +
      geom_rect(xmin=0,xmax=1,ymin=1,ymax=2,color="black",fill="blue",alpha=joints[3]) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,color="black",fill="blue",alpha=joints[4]) +
      annotate(geom="text",x=1.5,y=1.5,label=joints[1]) +
      annotate(geom="text",x=1.5,y=.5,label=joints[2]) +
      annotate(geom="text",x=.5,y=1.5,label=joints[3]) +
      annotate(geom="text",x=.5,y=.5,label=joints[4]) +
      scale_x_continuous(limits=c(0,2),breaks=c(.5,1.5),labels=c("Non-Effect","Effect")) +
      scale_y_continuous(limits=c(0,2),breaks=c(.5,1.5),labels=c("Non-Effect","Effect")) +
      labs(title=name,x=NULL,y=NULL) +
      theme_classic() +
      theme(aspect.ratio=1,plot.title=element_text(hjust=.5),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(angle=90,hjust=.5,vjust=-5),axis.text.x=element_text(vjust=6))
    return(plot)
  }

  w1 <- jointplot("sigma_w[1]","phi_w[1]",bquote(sigma[w]^{Area} ~ "x" ~ phi[w]^{Area}))
  w2 <- jointplot("sigma_w[2]","phi_w[2]",bquote(sigma[w]^{Altitude} ~ "x" ~ phi[w]^{Altitude}))
  e1 <- jointplot("sigma_e[1]","phi_e[1]",bquote(sigma[e]^{Area} ~ bold("*") ~ "x" ~ phi[e]^{Area} ~ bold("*")))
  e2 <- jointplot("sigma_e[2]","phi_e[2]",bquote(sigma[e]^{Altitude} ~ "x" ~ phi[e]^{Altitude}))
  b1 <- jointplot("sigma_b[1]","phi_b[1]",bquote(sigma[b]^{Distance} ~ "x" ~ phi[b]^{Distance}))
  b2 <- jointplot("sigma_b[2]","phi_b[2]",bquote(sigma[b]^{Altitude} ~ "x" ~ phi[b]^{Altitude}))
  d1 <- jointplot("sigma_d[1]","phi_d[1]",bquote(sigma[d]^{Distance} ~ bold("*") ~ "x" ~ phi[d]^{Distance} ~ bold("*")))
  d2 <- jointplot("sigma_d[2]","phi_d[2]",bquote(sigma[d]^{Altitude} ~ "x" ~ phi[d]^{Altitude}))

  axis_plot <- ggplot() +
    labs(x=bquote("Categorical (" ~ sigma ~ ")"),y=bquote("Quantitative (" ~ phi ~ ")")) +
    theme_classic() +
    theme(aspect.ratio=.6,line=element_blank())
  inner_plot <- w1 + w2 + e1 + e2 + b1 + b2 + d1 + d2 +
    plot_layout(ncol=4) &
    theme(text=element_text(size=8))
  joint_plot <- axis_plot + inset_element(inner_plot,left=0,bottom=0,right=1,top=1)

return(joint_plot)
}

make_statesplot <- function (file) {
  labels <- c("0" = "Aa","1" = "Ca","2" = "Pa","3" = "Cc","4" = "Ad","5" = "El","6" = "AaCa","7" = "AaPa","8" = "CaPa","9" = "AaCc","10" = "CaCc","11" = "PaCc","12" = "AaAd","13" = "CaAd","14" = "PaAd","15" = "CcAd","16" = "AaEl","17" = "CaEl","18" = "PaEl","19" = "CcEl","20" = "AdEl","21" = "AaCaPa","22" = "AaCaCc","23" = "AaPaCc","24" = "CaPaCc","25" = "AaCaAd","26" = "AaPaAd","27" = "CaPaAd","28" = "AaCcAd","29" = "CaCcAd","30" = "PaCcAd","31" = "AaCaEl","32" = "AaPaEl","33" = "CaPaEl","34" = "AaCcEl","35" = "CaCcEl","36" = "PaCcEl","37" = "AaAdEl","38" = "CaAdEl","39" = "PaAdEl","40" = "CcAdEl","41" = "AaCaPaCc","42" = "AaCaPaAd","43" = "AaCaCcAd","44" = "AaPaCcAd","45" = "CaPaCcAd","46" = "AaCaPaEl","47" = "AaCaCcEl","48" = "AaPaCcEl","49" = "CaPaCcEl","50" = "AaCaAdEl","51" = "AaPaAdEl","52" = "CaPaAdEl","53" = "AaCcAdEl","54" = "CaCcAdEl","55" = "PaCcAdEl","56" = "AaCaPaCcAd","57" = "AaCaPaCcEl","58" = "AaCaPaAdEl","59" = "AaCaCcAdEl","60" = "AaPaCcAdEl","61" = "CaPaCcAdEl","62" = "AaCaPaCcAdEl")
colors <- c("#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", "#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72")

states <- processAncStates(paste(lio_dir,"/",file,sep=""), state_labels=labels)

tree <- plotAncStatesMAP(t=states,
                   tree_layout="circular",
                   node_size=2,
                   node_color_as="state",
                   node_color=colors,
                   node_size_as=NULL) +
    ggplot2::theme(legend.direction="horizontal",
                   legend.position=c(.5,0),
                   legend.title=element_blank())

return(tree)
}

make_mdapeplot <- function () {
  all.errors <- data.frame(matrix(ncol=0,nrow=0))
  for (param in params) {
    param2 <- gsub("\\[|\\]", ".", param)
    true_colname <- paste(param2,".true",sep="")
    est_colname <- paste(param2,".mean",sep="")
    for (i in 1:nrow(regular_success_data)) {
      all.errors[i,"scenario"] <- paste(regular_success_data[i,"nregions"],regular_success_data[i,"treesize"],sep="")
      error <- (regular_success_data[i,est_colname]-regular_success_data[i,true_colname])/abs(regular_success_data[i,true_colname])*100
      all.errors[i,param2] <- error
    }
  }
  rep_str = c("xsmall"="X", "small"="S", "medium"="M", "large"="L")
  all.errors$scenario <- str_replace_all(all.errors$scenario, rep_str)

  all.errors_5X <- all.errors[which(all.errors$scenario=="5X"),]
  all.errors_5S <- all.errors[which(all.errors$scenario=="5S"),]
  all.errors_5M <- all.errors[which(all.errors$scenario=="5M"),]
  all.errors_5L <- all.errors[which(all.errors$scenario=="5L"),]
  all.errors_3M <- all.errors[which(all.errors$scenario=="3M"),]
  all.errors_7M <- all.errors[which(all.errors$scenario=="7M"),]

  rhos_5X <- as.vector(all.errors_5X[,2:5])
  all_rhos_5X <- c(rhos_5X$r_w,rhos_5X$r_e,rhos_5X$r_b,rhos_5X$r_d)
  sigmas_5X <- as.vector(all.errors[,6:13])
  all_sigmas_5X <- c(sigmas_5X$sigma_w.1.,sigmas_5X$sigma_w.2.,sigmas_5X$sigma_e.1.,sigmas_5X$sigma_e.2.,sigmas_5X$sigma_b.1.,sigmas_5X$sigma_b.2.,sigmas_5X$sigma_d.1.,sigmas_5X$sigma_d.2.)
  phis_5X <- as.vector(all.errors[,14:21])
  all_phis_5X <- c(phis_5X$phi_w.1.,phis_5X$phi_w.2.,phis_5X$phi_e.1.,phis_5X$phi_e.2.,phis_5X$phi_b.1.,phis_5X$phi_b.2.,phis_5X$phi_d.1.,phis_5X$phi_d.2.)

  rhos_5S <- as.vector(all.errors_5S[,2:5])
  all_rhos_5S <- c(rhos_5S$r_w,rhos_5S$r_e,rhos_5S$r_b,rhos_5S$r_d)
  sigmas_5S <- as.vector(all.errors_5S[,6:13])
  all_sigmas_5S <- c(sigmas_5S$sigma_w.1.,sigmas_5S$sigma_w.2.,sigmas_5S$sigma_e.1.,sigmas_5S$sigma_e.2.,sigmas_5S$sigma_b.1.,sigmas_5S$sigma_b.2.,sigmas_5S$sigma_d.1.,sigmas_5S$sigma_d.2.)
  phis_5S <- as.vector(all.errors_5S[,14:21])
  all_phis_5S <- c(phis_5S$phi_w.1.,phis_5S$phi_w.2.,phis_5S$phi_e.1.,phis_5S$phi_e.2.,phis_5S$phi_b.1.,phis_5S$phi_b.2.,phis_5S$phi_d.1.,phis_5S$phi_d.2.)

  rhos_5M <- as.vector(all.errors_5M[,2:5])
  all_rhos_5M <- c(rhos_5M$r_w,rhos_5M$r_e,rhos_5M$r_b,rhos_5M$r_d)
  sigmas_5M <- as.vector(all.errors_5M[,6:13])
  all_sigmas_5M <- c(sigmas_5M$sigma_w.1.,sigmas_5M$sigma_w.2.,sigmas_5M$sigma_e.1.,sigmas_5M$sigma_e.2.,sigmas_5M$sigma_b.1.,sigmas_5M$sigma_b.2.,sigmas_5M$sigma_d.1.,sigmas_5M$sigma_d.2.)
  phis_5M <- as.vector(all.errors_5M[,14:21])
  all_phis_5M <- c(phis_5M$phi_w.1.,phis_5M$phi_w.2.,phis_5M$phi_e.1.,phis_5M$phi_e.2.,phis_5M$phi_b.1.,phis_5M$phi_b.2.,phis_5M$phi_d.1.,phis_5M$phi_d.2.)

  rhos_5L <- as.vector(all.errors_5L[,2:5])
  all_rhos_5L <- c(rhos_5L$r_w,rhos_5L$r_e,rhos_5L$r_b,rhos_5L$r_d)
  sigmas_5L <- as.vector(all.errors_5L[,6:13])
  all_sigmas_5L <- c(sigmas_5L$sigma_w.1.,sigmas_5L$sigma_w.2.,sigmas_5L$sigma_e.1.,sigmas_5L$sigma_e.2.,sigmas_5L$sigma_b.1.,sigmas_5L$sigma_b.2.,sigmas_5L$sigma_d.1.,sigmas_5L$sigma_d.2.)
  phis_5L <- as.vector(all.errors_5L[,14:21])
  all_phis_5L <- c(phis_5L$phi_w.1.,phis_5L$phi_w.2.,phis_5L$phi_e.1.,phis_5L$phi_e.2.,phis_5L$phi_b.1.,phis_5L$phi_b.2.,phis_5L$phi_d.1.,phis_5L$phi_d.2.)

  rhos_3M <- as.vector(all.errors_3M[,2:5])
  all_rhos_3M <- c(rhos_3M$r_w,rhos_3M$r_e,rhos_3M$r_b,rhos_3M$r_d)
  sigmas_3M <- as.vector(all.errors_3M[,6:13])
  all_sigmas_3M <- c(sigmas_3M$sigma_w.1.,sigmas_3M$sigma_w.2.,sigmas_3M$sigma_e.1.,sigmas_3M$sigma_e.2.,sigmas_3M$sigma_b.1.,sigmas_3M$sigma_b.2.,sigmas_3M$sigma_d.1.,sigmas_3M$sigma_d.2.)
  phis_3M <- as.vector(all.errors_3M[,14:21])
  all_phis_3M <- c(phis_3M$phi_w.1.,phis_3M$phi_w.2.,phis_3M$phi_e.1.,phis_3M$phi_e.2.,phis_3M$phi_b.1.,phis_3M$phi_b.2.,phis_3M$phi_d.1.,phis_3M$phi_d.2.)

  rhos_7M <- as.vector(all.errors_7M[,2:5])
  all_rhos_7M <- c(rhos_7M$r_w,rhos_7M$r_e,rhos_7M$r_b,rhos_7M$r_d)
  sigmas_7M <- as.vector(all.errors_7M[,6:13])
  all_sigmas_7M <- c(sigmas_7M$sigma_w.1.,sigmas_7M$sigma_w.2.,sigmas_7M$sigma_e.1.,sigmas_7M$sigma_e.2.,sigmas_7M$sigma_b.1.,sigmas_7M$sigma_b.2.,sigmas_7M$sigma_d.1.,sigmas_7M$sigma_d.2.)
  phis_7M <- as.vector(all.errors_7M[,14:21])
  all_phis_7M <- c(phis_7M$phi_w.1.,phis_7M$phi_w.2.,phis_7M$phi_e.1.,phis_7M$phi_e.2.,phis_7M$phi_b.1.,phis_7M$phi_b.2.,phis_7M$phi_d.1.,phis_7M$phi_d.2.)

  order <- c("7M","5L","5M","5S","5X","3M")
  mapes <- data.frame(matrix(ncol=3,nrow=18))
  mapes[1:6,1] <- order
  mapes[7:12,1] <- order
  mapes[13:18,1] <- order
  mapes[1:6,2] <- rep("rho",6)
  mapes[7:12,2] <- rep("sigma",6)
  mapes[13:18,2] <- rep("phi",6)
  colnames(mapes) <- c("scenario","param","value")
  mapes[1,3] <- round(median(abs(all_rhos_7M)),2)
  mapes[2,3] <- round(median(abs(all_rhos_5L)),2)
  mapes[3,3] <- round(median(abs(all_rhos_5M)),2)
  mapes[4,3] <- round(median(abs(all_rhos_5S)),2)
  mapes[5,3] <- round(median(abs(all_rhos_5X)),2)
  mapes[6,3] <- round(median(abs(all_rhos_3M)),2)
  mapes[7,3] <- round(median(abs(all_sigmas_7M)),2)
  mapes[8,3] <- round(median(abs(all_sigmas_5L)),2)
  mapes[9,3] <- round(median(abs(all_sigmas_5M)),2)
  mapes[10,3] <- round(median(abs(all_sigmas_5S)),2)
  mapes[11,3] <- round(median(abs(all_sigmas_5X)),2)
  mapes[12,3] <- round(median(abs(all_sigmas_3M)),2)
  mapes[13,3] <- round(median(abs(all_phis_7M)),2)
  mapes[14,3] <- round(median(abs(all_phis_5L)),2)
  mapes[15,3] <- round(median(abs(all_phis_5M)),2)
  mapes[16,3] <- round(median(abs(all_phis_5S)),2)
  mapes[17,3] <- round(median(abs(all_phis_5X)),2)
  mapes[18,3] <- round(median(abs(all_phis_3M)),2)

  mdape_plot <- ggplot(mapes,aes(x=factor(scenario,level=order),y=value,group=param)) +
    geom_point(aes(color=param)) +
    geom_line(aes(color=param)) +
    scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20)) +
    scale_color_manual(name = "", values = c("rho"="red","phi"="green","sigma"="blue")) +
    labs(y="Median Absolute Percent Error (MDAPE)", x="Scenario") +
    theme_bw() +
    theme(aspect.ratio=.75,legend.position="top",plot.margin=margin(8,8,8,8))

  return(mdape_plot)
}

make_rateplots <- function() {
  order <- c("Aa","Ca","Pa","Cc","Ad","El")
  region_names <- eval(parse(text=paste("sf_data$","Region",sep="")))
  centroids <- gCentroid(sf_data, byid=TRUE)
  centroids <- as.data.frame(centroids)
  centroids$id <- region_names
  colnames(centroids) <- c("x","y","id")
  centroids$x <- lapply(centroids$x,function(z) round(z/100000))
  centroids$y <- lapply(centroids$y,function(z) round(z/100000))
  centroids <- centroids %>% slice(match(order,id))
  centroids$x <- as.numeric(centroids$x)
  centroids$y <- as.numeric(centroids$y)
  rownames(centroids) <- NULL
  centroids[6,1] <- centroids[6,1] - 8
  centroids[6,2] <- centroids[6,2] + 4
  centroids[5,1] <- centroids[5,1] + 2
  layout <- data.matrix(centroids[,1:2])

  d_weights <- matrix(NA,nrow=6,ncol=6)
  for (i in 1:6) {
    for (j in 1:6) {
      if (i == j) {
        d_weights[i,j] <- 0
      }
      else {
        colname <- paste("m_d.",i,"..",j,".",sep="")
        val <- mean(m_data[,colname])
        d_weights[i,j] <- val
      }
    }
  }
  colnames(d_weights) <- order
  d_weights <- d_weights

  b_weights <- matrix(NA,nrow=6,ncol=6)
  for (i in 1:6) {
    for (j in 1:6) {
      if (i == j) {
        b_weights[i,j] <- 0
      }
      else {
        colname <- paste("m_b.",i,"..",j,".",sep="")
        val <- mean(m_data[,colname])
        b_weights[i,j] <- val
      }
    }
  }
  colnames(b_weights) <- order
  b_weights <- b_weights

  e_sizes <- rep(NA,6)
  for (i in 1:6) {
    colname <- paste("m_e.",i,".",sep="")
    val <- mean(m_data[,colname])
    e_sizes[i] <- val
  }
  e_sizes <- sqrt(e_sizes) * 20

  w_rates <- rep(NA,6)
  w_sizes <- rep(NA,6)
  for (i in 1:6) {
    colname <- paste("m_w.",i,".",sep="")
    val <- mean(m_data[,colname])
    med <- median(m_data[,colname])
    w_sizes[i] <- val
    w_rates[i] <- med
  }
  w_rates <- data.frame(t(w_rates))
  colnames(w_rates) <- order
  w_sizes <- sqrt(w_sizes) * 20

  d_graph <- graph_from_adjacency_matrix(d_weights, mode="undirected", weighted=TRUE)
  b_graph <- graph_from_adjacency_matrix(b_weights, mode="undirected", weighted=TRUE)

  d_e_plot <- plot(d_graph,layout=layout,edge.width=E(d_graph)$weight,vertex.size=e_sizes,vertex.color="lightgrey",edge.color="black",vertex.label.color="black")
  de <- recordPlot()
  b_w_plot <- plot(b_graph,layout=layout,edge.width=E(b_graph)$weight,vertex.size=w_sizes,vertex.color="lightgrey",edge.color="black",vertex.label.color="black")
  bw <- recordPlot()
  return(list(de=de,bw=bw,w_rates=w_rates))
}

make_map <- function() {
  world <- ne_countries(scale = "medium", returnclass = "sf")
  sam <- subset(world, continent == "South America")
  box <- c(xmin=-100,xmax=-20,ymin=-50,ymax=20)
  sam <- st_crop(sam, box)
  region_names <- eval(parse(text=paste("sf_data$","Region",sep="")))
  plot <- ggplot(data=sam) +
    geom_sf() +
    coord_sf(crs = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs") +
    geom_polygon(data=fortify(sf_data), aes(x=long,y=lat,fill=id),alpha=.5) +
    scale_fill_manual(labels = region_names, values = c("red", "orange", "yellow", "green", "blue", "purple")) +
    labs(x="Longitude",y="Latitude",fill="Region") +
    theme_bw()
  return(plot)
}

make_onoffplot <- function () {
  sums <- rowSums(data_onoff)
  onoff_sums <- cbind(data_onoff, sums)
  sums <- rep(sums, 10)
  onoff <- data.frame(sums)
  total <- nrow(data_onoff)

  total_df <- data.frame(matrix(ncol=2,nrow=1))
  for (i in rjparams) {
    param <- paste("rj_",gsub("\\[|\\]", ".", i),sep="")
    count <- sum(data_onoff[,param])
    row <- data.frame(t(c(i,count)))
    total_df <- rbind(total_df,row)
  }
  total_df <- total_df[-(1),]
  colnames(total_df) <- c("param","percent")
  total_df$percent <- as.numeric(total_df$percent)/total

  stacked_df <- data.frame(matrix(ncol=3,nrow=1))
  for (i in 0:16) {
    onoff_sums_i <- onoff_sums[which(onoff_sums$sums==i),]
    for (j in rjparams) {
      param <- paste("rj_",gsub("\\[|\\]", ".", j),sep="")
      count <- sum(onoff_sums_i[,param])
      row <- data.frame(t(c(i,j,count)))
      stacked_df <- rbind(stacked_df,row)
    }
  }
  stacked_df <- stacked_df[-(1),]
  colnames(stacked_df) <- c("onoff","param","count")
  stacked_df$count <- as.numeric(stacked_df$count)

  levels <- c("phi_d[1]","sigma_d[1]","sigma_e[1]","phi_e[1]","sigma_e[2]","phi_b[1]","sigma_b[1]","sigma_b[2]","phi_e[2]","sigma_d[2]","phi_b[2]","phi_d[2]","sigma_w[2]","sigma_w[1]","phi_w[2]","phi_w[1]")
  labels <- c(
    bquote(phi[d]^{Distance} ~ bold(" * ")),
    bquote(sigma[d]^{Distance} ~ bold(" * ")),
    bquote(sigma[e]^{Area} ~ bold(" * ")),
    bquote(phi[e]^{Area} ~ bold(" * ")),
    bquote(sigma[e]^{Altitude}),
    bquote(phi[b]^{Distance}),
    bquote(sigma[b]^{Distance}),
    bquote(sigma[b]^{Altitude}),
    bquote(phi[e]^{Altitude}),
    bquote(sigma[d]^{Altitude}),
    bquote(phi[b]^{Altitude}),
    bquote(phi[d]^{Altitude}),
    bquote(sigma[w]^{Altitude}),
    bquote(sigma[w]^{Area}),
    bquote(phi[w]^{Altitude}),
    bquote(phi[w]^{Area}))

  palette1 <- brewer.pal(8,"Dark2")
  palette2 <- brewer.pal(8,"Set2")
  palette3 <- c(rbind(palette1,palette2))

  total_plot <- ggplot(total_df,aes(x=factor(param,levels=levels),y=percent,fill=factor(param,levels=levels))) +
    geom_bar(stat="identity") +
    scale_fill_manual(labels=labels,values=palette3) +
    scale_y_continuous(limits=c(0,1),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0)) +
    labs(fill="Parameter",x="Parameter",y="Frequency") +
    theme_bw() +
    theme(aspect.ratio=.4,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(8,8,4,8))

  onoff_plot <- ggplot(onoff, aes(sums)) +
    geom_histogram(aes(y=..density..), binwidth=1, boundary=.5, color="black", size=.1, fill="white") +
    geom_density(aes(linetype="Observed"), adjust=4, linewidth=1) +
    geom_density(aes(rbinom(20010,16,.5), linetype="Prior"), adjust=4, linewidth=1) +
    scale_linetype_manual(values=c("Observed"="solid","Prior"="dashed")) +
    scale_x_continuous(breaks=seq(0,16,1), limits=c(-.5,16.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,.3),breaks=c(0.0,0.1,0.2,0.3),labels=c("0.00","0.10","0.20","0.30"),expand=c(0,0)) +
    labs(y="Density", x=NULL, linetype="Distribution") +
    theme_bw() +
    theme(aspect.ratio=.3,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(4,8,4,8))

  stacked_plot <- ggplot(stacked_df,aes(x=onoff,y=count,fill=factor(param,levels=levels))) +
    geom_bar(stat="identity",position="fill",width=1,color="black",size=.1) +
    scale_x_discrete(limits=factor(seq(0,16,1)),expand=c(0,0)) +
    scale_y_continuous(breaks=seq(0,1,.25),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0)) +
    scale_fill_manual(labels=labels,values=palette3) +
    labs(fill="Parameter",x="Number of 'ON' Parameters",y="Representation By Bin") +
    theme_bw() +
    theme(aspect.ratio=.75,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none",plot.margin=margin(4,8,8,8))

  big_onoff_plot <- total_plot + onoff_plot + stacked_plot + plot_layout(ncol=1,guides="collect")

  return(big_onoff_plot)
}

##################################
########## DATA LOADING ##########
##################################

# File Directories
in_dir <- "../../output/processed_data"
lio_dir <- "../../output/emp"
out_dir <- "../../output/plots"

# Miscellaneous Data
geo_variances_data <- read.csv(file=paste(in_dir, "/", "geo_variances.csv", sep=""),header=TRUE)
lio_variances_data <- read.csv(file=paste(in_dir, "/", "lio_variances.csv", sep=""),header=TRUE)
sf_data <- raster::shapefile("../../data/emp/shapefiles/sam_regions_final.shp")

# Parameters
params <- c("r_w","r_e","r_b","r_d","sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]")
param_names <- c("r_w","r_e","r_b","r_d","sigma_w[Area]","sigma_w[Altitude]","sigma_e[Area]","sigma_e[Altitude]","sigma_b[Distance]","sigma_b[Altitude]","sigma_d[Distance]","sigma_d[Altitude]","phi_w[Area]","phi_w[Altitude]","phi_e[Area]","phi_e[Altitude]","phi_b[Distance]","phi_b[Altitude]","phi_d[Distance]","phi_d[Altitude]")
rjparams <- c("sigma_w[1]","sigma_w[2]","sigma_e[1]","sigma_e[2]","sigma_b[1]","sigma_b[2]","sigma_d[1]","sigma_d[2]","phi_w[1]","phi_w[2]","phi_e[1]","phi_e[2]","phi_b[1]","phi_b[2]","phi_d[1]","phi_d[2]")

#############################################
########## SIM COMMANDS AND SAVING ##########
#############################################

# Sim Data
regular_sim_data <- read.csv(file=paste(in_dir, "/", "regular_sim.csv", sep=""),header=TRUE)
rj_sim_data <- read.csv(file=paste(in_dir, "/", "rj_sim.csv", sep=""),header=TRUE)
regular_success_data <- subset(regular_sim_data, success==TRUE)
m3 <- subset(regular_sim_data, nregions==3 & treesize=="medium" & success==TRUE)
x5 <- subset(regular_sim_data, nregions==5 & treesize=="xsmall" & success==TRUE)
s5 <- subset(regular_sim_data, nregions==5 & treesize=="small" & success==TRUE)
m5 <- subset(regular_sim_data, nregions==5 & treesize=="medium" & success==TRUE)
l5 <- subset(regular_sim_data, nregions==5 & treesize=="large" & success==TRUE)
m7 <- subset(regular_sim_data, nregions==7 & treesize=="medium" & success==TRUE)
m5rj <- subset(rj_sim_data, nregions==5 & treesize=="medium" & success==TRUE)

# Map
map <- make_map()
pdf(file=paste(out_dir, "/map.pdf", sep=""))
map
dev.off()

# Pairwise Variance Plot
varplot <- make_varplot()
pdf(file=paste(out_dir, "/variances.pdf", sep=""))
varplot
dev.off()

# Coverage Plots (plot & table)
covplots <- make_covplots()
pdf(file=paste(out_dir, "/coverageplot.pdf", sep=""))
covplots$plot
dev.off()
gtsave(covplots$table,paste(out_dir, "/coverages.pdf", sep=""))

# RJ Sim Plot
rjplot <- make_rjplot()
pdf(file=paste(out_dir, "/rjprobplot.pdf", sep=""))
rjplot
dev.off()

# MDAPE Plot
mdapeplot <- make_mdapeplot()
pdf(file=paste(out_dir, "/mdape.pdf", sep=""))
mdapeplot
dev.off()

############################################
########## HL COMMANDS AND SAVING ##########
############################################

# Emperical Data
data_nonrj <- read.csv(file=paste(lio_dir,"/",51,".model.log",sep=""),sep ="\t",header=TRUE)
data_rj <- read.csv(file=paste(lio_dir,"/",52,".model.log",sep=""),sep="\t",header=TRUE)
m_data <- read.csv(file=paste("~/projects/bg_liolaemidae/emp_output/51.geo_features.log"),header=TRUE,sep="\t")
data_onoff <- data_rj[,161:176]
sums <- rowSums(data_onoff)
prop <- sum(sums>4)/length(sums)

# Posterior Plots (r, w, e, b, d)
postplots <- make_postplots()
pdf(file=paste(out_dir, "/hl_w.pdf", sep=""))
postplots$w
dev.off()
pdf(file=paste(out_dir, "/hl_e.pdf", sep=""))
postplots$e
dev.off()
pdf(file=paste(out_dir, "/hl_b.pdf", sep=""))
postplots$b
dev.off()
pdf(file=paste(out_dir, "/hl_d.pdf", sep=""))
postplots$d
dev.off()
pdf(file=paste(out_dir, "/hl_baserates.pdf", sep=""))
postplots$r
dev.off()

# Joint Plot
jointplot <- make_jointplot()
pdf(file=paste(out_dir, "/hl_jointprobs.pdf", sep=""))
jointplot
dev.off()

# Ancestral States
statesplot <- make_statesplot("51.ase.tre")
ggsave(paste(out_dir,"/hl_states.pdf",sep=""),width=9,height=12)

# Rate Plots
rateplots <- make_rateplots()
pdf(file=paste(out_dir, "/hl_derates.pdf", sep=""))
rateplots$de
dev.off()
pdf(file=paste(out_dir, "/hl_bwrates.pdf", sep=""))
rateplots$bw
dev.off()
rateplots$w_rates

# On/Off Plot
onoffplot <- make_onoffplot()
pdf(file=paste(out_dir, "/hl_onoff.pdf", sep=""))
onoffplot
dev.off()

############################################
########## AN COMMANDS AND SAVING ##########
############################################

# Emperical Data
data_nonrj <- read.csv(file=paste(lio_dir,"/",53,".model.log",sep=""),sep ="\t",header=TRUE)
data_rj <- read.csv(file=paste(lio_dir,"/",54,".model.log",sep=""),sep="\t",header=TRUE)
m_data <- read.csv(file=paste("~/projects/bg_liolaemidae/emp_output/53.geo_features.log"),header=TRUE,sep="\t")
data_onoff <- data_rj[,161:176]

# Posterior Plots (r, w, e, b, d)
postplots <- make_postplots()
pdf(file=paste(out_dir, "/an_w.pdf", sep=""))
postplots$w
dev.off()
pdf(file=paste(out_dir, "/an_e.pdf", sep=""))
postplots$e
dev.off()
pdf(file=paste(out_dir, "/an_b.pdf", sep=""))
postplots$b
dev.off()
pdf(file=paste(out_dir, "/an_d.pdf", sep=""))
postplots$d
dev.off()
pdf(file=paste(out_dir, "/an_baserates.pdf", sep=""))
postplots$r
dev.off()

# Joint Plot
jointplot <- make_jointplot()
pdf(file=paste(out_dir, "/an_jointprobs.pdf", sep=""))
jointplot
dev.off()

# Ancestral States
statesplot <- make_statesplot("53.ase.tre")
ggsave(paste(out_dir,"/an_states.pdf",sep=""),width=9,height=12)

# Rate Plots
rateplots <- make_rateplots()
pdf(file=paste(out_dir, "/an_derates.pdf", sep=""))
rateplots$de
dev.off()
pdf(file=paste(out_dir, "/an_bwrates.pdf", sep=""))
rateplots$bw
dev.off()
rateplots$w_rates

# On/Off Plot
onoffplot <- make_onoffplot()
pdf(file=paste(out_dir, "/an_onoff.pdf", sep=""))
onoffplot
dev.off()

######################################
########## REDUCED FEATURES ##########
######################################

# Emperical Data
data_nonrj <- read.csv(file=paste(lio_dir,"/",55,".model.log",sep=""),sep ="\t",header=TRUE)
data_rj <- read.csv(file=paste(lio_dir,"/",56,".model.log",sep=""),sep="\t",header=TRUE)
m_data <- read.csv(file=paste("~/projects/bg_liolaemidae/emp_output/55.geo_features.log"),header=TRUE,sep="\t")
data_onoff <- data_rj[,161:176]

# Joint Plot
jointplot <- make_jointplot()
pdf(file=paste(out_dir, "/rf_jointprobs.pdf", sep=""))
jointplot
dev.off()

# On/Off Plot
onoffplot <- make_onoffplot()
pdf(file=paste(out_dir, "/rf_onoff.pdf", sep=""))
onoffplot
dev.off()
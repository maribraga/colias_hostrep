#'---
#'title: "Colias diversification"
#'author: "Mariana P Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#+ setup, include = FALSE
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

#'-------------
#'
#' Script for analyses performed in ... 
#' *...*
#'

#' ## Set up
#' For this script we'll need  

#+ packages
library(evolnets)
library(MCMCpack)
library(coda)
library(kdensity)
library(igraph)
library(ape)
library(treeio)
library(ggtree)
library(patchwork)
library(bipartite)
library(tidyverse)

  
#' ## 
#'

#' ** **
#' 
# read files
path_data <- "./diversification/data/"


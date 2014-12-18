library(ggplot2)
library(reshape2)
library(RColorBrewer)
#install.packages('devtools')
#require(devtools)
#install_github('likert','jbryer')
library('likert')
require(devtools)

#' Calculate the likert
#' 
#' @param cluster_result Clustering result eg: [1, 2, 1, 3 ,3]
#' @param sample_stage Stage in which the sample is eg: [3, 3, 1, 2 ,2]
#' @param cluster_prefix Prefix to prefend to cluster data
#' @param sample_prefix Prefix to prefend to stage data
#' @return The sum of \code{x} and \code{y}.
likertToClus <- function(cluster_result, sample_stage, cluster_prefix='K', sample_prefix='S'){
  # check different value
  cluster <- unique(cluster_result)
  stage <- unique(sample_stage)
  
  # create label ['1', '2'] + prefix='k' -> ['K1', 'K2'] 
  cluster_label <- paste(cluster_prefix, cluster, sep='')
  stage_label <- paste(sample_prefix, stage, sep='')
  
  # create factor based on clustering and stage data
  factor_k <- factor(cluster_result, labels = cluster_label)
  factor_s <- factor(sample_stage, labels = stage_label)
  
  # data frame c1->cluster, c2->stage
  df <- data.frame(factor_k, factor_s)
  col <- df[,2, drop=F]
  # calc likert  
  result <- likert(col, grouping = df$factor_k)
  result
}

test_likertToClus <- function() {
  # sample data
  K <- 4
  S <- 6
  sample <- 250
  example_label <- sample(1:K, size = sample, replace = T)
  example_stage <- sample(1:S, size = sample, replace = T)

  # exec
  tab <- likertToClus(example_label, example_stage)
  plot(tab, centered = F, include.histogram = TRUE, col = terrain.colors(S))
}

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
  cluster <- sort(unique(cluster_result), decreasing = T)
  stage <- sort(unique(sample_stage), na.last = T)
  
  # create label ['1', '2'] + prefix='k' -> ['K1', 'K2'] 
  cluster_label <- paste(cluster_prefix, cluster, sep='')
  stage_label <- paste(sample_prefix, stage, sep='')
  
  # create factor based on clustering and stage data
  factor_k <- factor(cluster_result, labels = cluster_label, exclude = NULL)
  factor_s <- factor(sample_stage, labels = stage_label, exclude = NULL)
  
  # data frame c1->cluster, c2->stage
  df <- data.frame(factor_k, factor_s)
  col <- df[,2, drop=F]
  # calc likert  
  result <- likert(col, grouping = df$factor_k)
  result
}

test_likertToClus <- function(K = 4, S = 6, sample = 250) {
  # sample data
  example_label <- sample(1:K, size = sample, replace = T)
  example_stage <- sample(1:S, size = sample, replace = T)

  # exec
  tab <- likertToClus(example_label, example_stage)
  plot(tab, centered = F, include.histogram = TRUE, col = terrain.colors(S))
}


# TODO: il blocco successivo va in una funzione


sample_NBS_cluster <- read.csv('/home/dex/Dropbox/pyTSA/TCGA-CRC likert/samples-NBS-mapping.txt', header = TRUE, sep = ' ')
sample_NBS_data <- read.csv('/home/dex/Dropbox/pyTSA/TCGA-CRC likert/tcga_infoClin_paper.txt', 
                            header = TRUE, 
                            sep = '\t', 
                            row.names = 'patient')

tumor_stage <-  subset(sample_NBS_data, select = 'tumor_stage')
sample_NBS_cluster$old_name <- sample_NBS_cluster$Sample.ID
sample_NBS_cluster$Sample.ID <- substr(sample_NBS_cluster$Sample.ID, 1, 12)
sample_NBS_cluster <- sample_NBS_cluster[!duplicated(sample_NBS_cluster$Sample.ID), ]
tumor_cluster = data.frame(sample_NBS_cluster$Cluster, sample_NBS_cluster$old_name, row.names=sample_NBS_cluster$Sample.ID)
sample_stage_cluster <- merge(tumor_stage, tumor_cluster, by = 0, all=T)

#levels(sample_stage_cluster$tumor_stage)
levels(sample_stage_cluster$sample_NBS_cluster.Cluster) <- sort(unique(sample_stage_cluster$sample_NBS_cluster.Cluster))
#levels(sample_stage_cluster$sample_NBS_cluster.Cluster)

result <- likertToClus(sample_stage_cluster$sample_NBS_cluster.Cluster, sample_stage_cluster$tumor_stage, sample_prefix = '', cluster_prefix ='')
plot(result, centered=F, include.histogram = TRUE, col = brewer.pal(length(unique(sample_stage_cluster$tumor_stage)), "YlOrRd"))

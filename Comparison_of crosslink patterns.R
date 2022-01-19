######################################################################
###
### Comparison of crosslink patterns in different iCLIP experiments
### 
###
######################################################################

### ===============================================================
### Structure of this script
### ---------------------------------------------------------------
###
### 1) Normalisation and spline smoothing of crosslink patterns
### 2) smoothed crosslink patterns around binding sites
### 3) crosslink patterns in the start of 3'UTRs
###
###

### ===============================================================
### 0) Libraries & Input
### ===============================================================

# Libraries
library(tidyverse)
library(GenomicRanges)
library(ggrastr)
library(Cairo)
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

# Input
BS # a Granges object of binding sites definde from endogenous PURA iCLIP
CL_path_experimentwise # a list of paths to the bigwig files containing the crosslink events per experiment (here split in plus and minus strand)
anno_full # gff3 file of Genome annotation from GENCODE imported as Granges object

### ===============================================
### 1) Normalisation and spline smoothing of crosslink patterns
### ===============================================
# Note: this approach was adapted from Heyl & Backofen 2021 - "StoatyDive: Evaluation and classification of peak profiles for sequencing data"

# min max norm
min_max_noramlize <- function(x){
  x_new <- (x-min(x)) / (max(x)-min(x)) 
  return(x_new)
}


# smooth
smoothing <- function(y, lambda, dim){
  data_points <- length(y)
  
  # Test if we have just a vector of constant values.
  which(y != max(y))
  if ( length(which(y != max(y))) != 0 ){
    x <- c(1:data_points)
    y.smooth <- smooth.spline(x, y, spar=lambda)
    x_new <- seq(1, data_points, data_points/dim)
    y_new <- predict(y.smooth, x_new)$y
    y_new <- y_new/max(y_new)
    return(y_new)
  } 
  
  x_new <- seq(1, data_points, data_points/dim)
  y_new <- rep( max(y), length(x_new) )
  return(y_new)
}


clean_smooth <- function(data_smoothed) {
  
  data_smoothed <- t(data_smoothed)
  # get rid of negative values
  data_smoothed[which(data_smoothed < 0)] = 0.0
  
  # get rid of diract impulses bigger than max (> 1.0)
  data_smoothed[which(data_smoothed > 1.0)] = 1.0 
  
  return(data_smoothed)
}


### ====================================
### 2) smoothed crosslink patterns around binding sites
### ====================================
window_around_BS # the size by which BS are enlarged to both sites to obtain the window of interest

# import crosslinks
crosslink_rles <- map(CL_path_experimentwise, ~import.bw(.x, as="RleList"))
cl_window <- BS + window_around_BS
cl_window <- c(map(crosslink_rles[c(1,3,5)], ~.x[window[strand(window)=="+"]] %>% as.matrix() ),
               map(crosslink_rles[c(2,4,6)], ~.x[window[strand(window)=="-"]] %>% as.matrix() %>% .[,ncol(.):1] ))

# normalise
cl_window_normalized <- lapply(cl_window, function(x) apply(x, 1, min_max_noramlize) %>% t(.))


# Set NaNs to zero
cl_window_normalized[[1]][which(is.na(cl_window_normalized[[1]]))] = 0
cl_window_normalized[[2]][which(is.na(cl_window_normalized[[2]]))] = 0
cl_window_normalized[[3]][which(is.na(cl_window_normalized[[3]]))] = 0

# smooth
cl_window_norm_smooth <- lapply(cl_window_normalized,  apply, 1, function(x) smoothing(y=x, lambda = 0.25, dim =100))
cl_window_norm_smooth <- lapply(cl_window_norm_smooth, function(x) clean_smooth(data_smoothed = x))  


# cluster within heatmap
cluster <- kmeans(cl_window_norm_smooth, 4)

cl_window_norm_smooth_cluster <- data.frame(cl_window_norm) %>% 
  mutate(., exp = c(rep("endo", 1000), rep("oe", 1000), rep("flag", 1000)),
         cluster = as.character(cluster$cluster),
         BS_region = rep(window$BS_region, times =3),
         BS_idx = rep(times =3, 1:1000)) %>% 
  arrange(cluster)

# make heatmap
Heatmap(as.matrix(cl_window_norm_smooth_cluster[,1:94], ncol = 94), cluster_columns = F, cluster_rows = F,
        right_annotation = rowAnnotation( experiment =  cl_window_norm_smooth_cluster$exp , 
                                          BS_region =  cl_window_norm_smooth_cluster$BS_region,
                                          col = list(experiment = c("endo" = "red", "oe" = "green", "flag" = "blue"),
                                                     BS_region = c("CDS" = "red", "3UTR" = "blue", "5UTR" = "darkgrey", "intron" = "grey", "non_cod" = "white"))), 
        row_split = cl_window_norm_smooth_cluster$cluster)



### ============================================
### 3) crosslink patterns in the start of 3'UTRs
### ============================================

# window at the beginning of the 3'UTR
# select UTRs
three_UTRs <- anno_full[anno_full$type=="three_prime_UTR"] 
three_UTRs <- three_UTRs[!duplicated(three_UTRs$gene_id)]
w = 300

window <- three_UTRs
end(window[strand(window)=="+"]) <- start(window[strand(window)=="+"])+300
start(window[strand(window)=="-"]) <- end(window[strand(window)=="-"])-300


cl_window <- c(map(crosslink_rles[c(1,3,5)], ~.x[window[strand(window)=="+"]] %>% as.matrix() ),
               map(crosslink_rles[c(2,4,6)], ~.x[window[strand(window)=="-"]] %>% as.matrix() %>% .[,ncol(.):1] ))

# select 3'UTRs with medium coverage
cl_window_coverage <- map(cl_window, ~ rowSums(.x))
map(cl_window_coverage, ~hist(log(.x)))

cov_p <- map_dfc(cl_window_coverage[1:3], ~.x) %>% mutate(medium_cov =
                                                            case_when(`...1` > 100 &
                                                                        `...2` > 100 &
                                                                        `...3` > 100 &
                                                                        `...1` < 10^6 &
                                                                        `...2` < 10^6 &
                                                                        `...3` < 10^6 
                                                                      ~ T,
                                                                      T~F))
cov_m <- map_dfc(cl_window_coverage[4:6], ~.x) %>% mutate(medium_cov =
                                                            case_when(`...1` > 100 &
                                                                        `...2` > 100 &
                                                                        `...3` > 100 &
                                                                        `...1` < 10^6 &
                                                                        `...2` < 10^6 &
                                                                        `...3` < 10^6 
                                                                      ~ T,
                                                                      T~F))


cl_window <- c(map(cl_window[1:3], ~.x[cov_p$medium_cov,] ),
               map(cl_window[4:6], ~.x[cov_m$medium_cov,] ))


window_p <- window[strand(window)=="+"]
window_m <- window[strand(window)=="-"]
window <- c(window_p[cov_p$medium_cov], window_m[cov_m$medium_cov])


cl_window <- list(endo = rbind(cl_window[[2]], cl_window[[5]]),
                  oe  = rbind(cl_window[[1]], cl_window[[4]]),
                  flag = rbind(cl_window[[3]], cl_window[[6]]))


# normalise
cl_window_normalized <- lapply(cl_window, function(x) apply(x, 1, min_max_noramlize) %>% t(.))


# Set NaNs to zero
cl_window_normalized[[1]][which(is.na(cl_window_normalized[[1]]))] = 0
cl_window_normalized[[2]][which(is.na(cl_window_normalized[[2]]))] = 0
cl_window_normalized[[3]][which(is.na(cl_window_normalized[[3]]))] = 0

# smoothing
cl_window_norm_smooth <- lapply(cl_window_normalized,  apply, 1, function(x) smoothing(y=x, lambda = 0.5, dim =500))
cl_window_norm_smooth <- lapply(cl_window_norm_smooth, function(x) clean_smooth(data_smoothed = x))  


# Heatmap of all three
cl_window_norm_smooth <- rbind(cl_window_norm_smooth$endo, cl_window_norm_smooth$oe, cl_window_norm_smooth$flag)
h1 <- Heatmap(cl_window_norm_smooth[[1]], cluster_columns = F, use_raster = T, raster_device = "png", raster_quality = 10) 
h2 <- Heatmap(cl_window_norm_smooth[[2]], cluster_columns = F, use_raster = T, raster_device = "png", raster_quality = 10) 
h3 <- Heatmap(cl_window_norm_smooth[[3]], cluster_columns = F, use_raster = T, raster_device = "png", raster_quality = 10)  

h1+h2+h3
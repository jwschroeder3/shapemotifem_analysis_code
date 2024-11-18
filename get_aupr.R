#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)

MOTIF_LEN = 10

# written by chatgpt, edited by me
extract_min_distances <- function(distances) {
  m <- dim(distances)[3]
  r <- dim(distances)[2]
  
  #mdist <- numeric(r)

  mdist_l = apply(distances, FUN=min, MAR=2:3)
  perc_dist = apply(mdist_l, FUN=percent_rank, MAR=2)
  mdist = apply(perc_dist, FUN=min, MAR=1)
  
  #for (j in 1:r) {
  #  min_value <- min(distances[ , j , ])
  #  mdist[j] <- min_value
  #}
  
  return(mdist)
}

# written by chatgpt, edited by me
calculate_hits <- function(mdist) {
  min_value <- min(mdist)
  max_value <- max(mdist)
  
  thresholds <- seq(min_value, max_value, length.out = 100)
  hits <- matrix(0, nrow = length(thresholds), ncol = length(mdist))
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    hits[i, ] <- as.integer(mdist <= threshold)
  }
  
  return(hits)
}

# written by chatgpt, edited by me
calculate_precision_recall <- function(hits, scores) {
  if (ncol(hits) != length(scores)) {
    stop("Number of columns in hits must be the same as the length of scores")
  }
  
  thresholds <- nrow(hits)
  precision <- numeric(thresholds)
  recall <- numeric(thresholds)
  
  for (i in 1:thresholds) {
    tp <- sum(hits[i, ] == 1 & scores == 1)  # True positives
    fp <- sum(hits[i, ] == 1 & scores == 0)  # False positives
    fn <- sum(hits[i, ] == 0 & scores == 1)  # False negatives
    
    precision[i] <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
    recall[i] <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
  }
  
  return(data.frame(threshold_index = 1:thresholds, precision = precision, recall = recall))
}

# written by chatgpt, edited by me
te_aupr <- function(precision_recall) {
  precision <- precision_recall$precision
  recall <- precision_recall$recall

  # Sort by recall for integration
  sorted_indices <- order(recall)
  precision <- precision[sorted_indices]
  recall <- recall[sorted_indices]

  # Calculate the area under the precision-recall curve using the trapezoidal rule
  aupr <- 0
  for (i in 2:length(recall)) {
    aupr <- aupr + (recall[i] - recall[i - 1]) * (precision[i] + precision[i - 1]) / 2
  }
  
  return(aupr)
}

# Define the options
option_list = list(
  make_option(c("-s", "--std_data"), type = "character", default = NULL, 
              help = "Basename of the RData file containing scores and standardized target test sequences and standardized motifs", metavar = "file"),
  make_option(c("--distances"), type = "character", default = NULL, 
              help = "RData file containing manhattan distances", metavar = "file"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL, 
              help = "Characters to prepend to output files"),
  make_option(c("-d", "--data_direc"), type = "character", default = NULL, 
              help = "Directory to search for input files", metavar = "directory")
)

# Parse the command-line arguments
opt_parser = OptionParser(option_list = option_list)
options = parse_args(opt_parser)

# Check for missing arguments
if (is.null(options$std_data)) {
  stop("Error: --std_data is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$distances)) {
  stop("Error: --distances is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$prefix)) {
  stop("Error: --prefix is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$data_direc)) {
  stop("Error: --data_direc is required. Use --help for usage information.", call. = FALSE)
}

prefix = options$prefix

# Now you can use options$data_direc and options$np_file in your script:
setwd(options$data_direc)
load(options$std_data)
load(options$distances)
score_vec = converted$scores$score

min_dists = extract_min_distances(distances) 
hits = calculate_hits(min_dists)
prec_recall = calculate_precision_recall(hits, score_vec)
aupr = te_aupr(prec_recall) 

outname = paste0(prefix, "_aupr.txt")
con = file(outname)
writeLines(format(aupr), con)
close(con)
print(paste0("AUPR: ", aupr))
print(paste0("AUPR written to ", outname))


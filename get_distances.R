#!/usr/bin/env Rscript

library(ShapeMotifEM)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(optparse)

MOTIF_LEN = 10

# this function written by chatGPT then modified by me.
# query is the motifs, reference is the sequences
calculate_manhattan_distance <- function(query, reference) {
  # s is shape number
  s <- dim(query)[1]
  # n is the motif length
  n <- dim(query)[2]
  # m is the number of motifs
  m <- dim(query)[3]
  # l is the reference sequence lengths
  l <- dim(reference)[2]
  # r is the number of reference sequences
  r <- dim(reference)[3]

  if (s != dim(reference)[1]) {
    stop("s dimension must be the same for both query and reference")
  }

  distances <- array(NA, dim = c(l - n + 1, r, m))

  # for reference sequence j
  for (j in 1:r) {
    # for ref seq start position i
    for (i in 1:(l - n + 1)) {
      # for motif k
      for (k in 1:m) {
        # grab correct length reference segment from sequence j, grab all shapes
        segment <- reference[, i:(i + n - 1), j]
        # get motif k's manhattan distance
        distance <- sum(abs(query[,,k] - segment))
        # place motif k's manhattan distance into distances array for position i, ref j
        distances[i, j, k] <- distance
      }
    }
  }
  
  return(distances)
}

# Define the options
option_list = list(
  make_option(c("-s", "--std_data"), type = "character", default = NULL, 
              help = "Basename of the RData file containing scores and standardized target test sequences and standardized motifs", metavar = "file"),
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

# names of items in the list "converted" are "scores", "z_target_shapes", and "z_motifs"
distances = calculate_manhattan_distance(converted$z_motifs, converted$z_target_shapes)
save(distances, file=paste0(prefix, "_distances.RData"))

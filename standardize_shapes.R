#!/usr/bin/env Rscript

library(ShapeMotifEM)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(optparse)

MOTIF_LEN = 10

# Define the options
option_list = list(
  make_option(c("-m", "--motifs_file"), type = "character", default = NULL, 
              help = "Basename of the motifs RData file", metavar = "file"),
  make_option(c("-f", "--fa_file"), type = "character", default = NULL, 
              help = "Basename of the fasta file", metavar = "file"),
  make_option(c("-s", "--score_file"), type = "character", default = NULL, 
              help = "Basename of the file containing peak/non-peak assignment for each sequence in --fa_file", metavar = "file"),
  make_option(c("-d", "--data_direc"), type = "character", default = NULL, 
              help = "Directory to search for input files", metavar = "directory")
)

# Parse the command-line arguments
opt_parser = OptionParser(option_list = option_list)
options = parse_args(opt_parser)

# Check for missing arguments
if (is.null(options$motifs_file)) {
  stop("Error: --motifs_file is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$score_file)) {
  stop("Error: --score_file is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$fa_file)) {
  stop("Error: --fa_file is required. Use --help for usage information.", call. = FALSE)
}

if (is.null(options$data_direc)) {
  stop("Error: --data_direc is required. Use --help for usage information.", call. = FALSE)
}

base = basename(options$fa_file)
prefix = tools::file_path_sans_ext(base)

# Now you can use options$data_direc and options$np_file in your script:
setwd(options$data_direc)
load(file=options$motifs_file)

###################################
###################################
## I need to re-organize the motifs to be a matrix of shape (nShapes, motifLength, nMotifs)
## They are currently nested list. Top level of list is length=nShapes, next
## level of lis is length=nMotifs
###################################
###################################
motif_mat = array(
    0.0,
    dim=c(
        length(area_dfm_set),
        nrow(area_dfm_set[[1]][[1]]),
        length(area_dfm_set[[1]])
    )
)
z_motif_mat = array(
    0.0,
    dim=c(
        length(area_dfm_set),
        nrow(area_dfm_set[[1]][[1]]),
        length(area_dfm_set[[1]])
    )
)
for (i in 1:length(area_dfm_set)) {
    shape_vals = area_dfm_set[[i]]
    for (j in 1:length(shape_vals)) {
        # for the i-th shape and j-th motif,
        # get the mean value for each position
        motif_mat[i, , j] = shape_vals[[j]]$meanv
    }
}

# fa_input_result is array of shape (nShapes, seqLen, numSeqs)
fa_input_result = fasta_input(options$fa_file, shapeind = "all")
z_scores = array(0.0, dim=dim(fa_input_result))
for (i in 1:dim(fa_input_result)[1]) {
    mean_i = mean(fa_input_result[i, , ])
    sd_i = sd(fa_input_result[i, , ])
    # get z-score for target sequences
    z_scores[i, , ] = (fa_input_result[i, , ] - mean_i) / sd_i
    # get z-score for motifs
    z_motif_mat[i, , ] = (motif_mat[i, , ] - mean_i) / sd_i
}
scores = read_tsv(options$score_file)

converted = list(
    "scores" = scores,
    "z_target_shapes" = z_scores,
    "z_motifs" = z_motif_mat
)

fname = paste0(prefix, "_standardized_data.RData")
print(paste0("Writing standardized data to ", fname))
save(converted, file=fname)


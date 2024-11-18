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
#  make_option(c("-n", "--np_file"), type = "character", default = NULL, 
#              help = "Basename of the narrowpeak file", metavar = "file"),
  make_option(c("-f", "--fa_file"), type = "character", default = NULL, 
              help = "Basename of the fasta file", metavar = "file"),
  make_option(c("-d", "--data_direc"), type = "character", default = NULL, 
              help = "Directory to search for input files", metavar = "directory")
)

# Parse the command-line arguments
opt_parser = OptionParser(option_list = option_list)
options = parse_args(opt_parser)

# Check for missing arguments
#if (is.null(options$np_file)) {
#  stop("Error: --np_file is required. Use --help for usage information.", call. = FALSE)
#}

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

#extraCols_narrowPeak = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
#peak_gr = import(options$np_file, format="BED", extraCols=extraCols_narrowPeak)
#
#sinfo = seqinfo(Hsapiens)
#seqlevels(peak_gr) = seqlevels(Hsapiens)
#seqlengths(peak_gr) = seqlengths(Hsapiens)
#seqinfo(peak_gr) = sinfo
#
## This uses BSgenome.Hsapiens.UCSC.hg19 on the backend, so won't work
## for many applications. Maybe I can instead just get straight to the result, which
## is just an array of shape (shape_num, peak_length, peak_num)
#gr_input_result = gr_input(
#  peak_gr,
#  # how much of center of peak to grab
#  peakLength = 100,
#  # how many peaks to grab after sorting
#  total_peak_number = 100,
#  # output fasta names
#  fasta_name = "fasta.fa",
#  shapeind = "all",
#  sort_by = "signalValue",
#  sort_descend = TRUE
#)
#
#save(gr_input_result, file="gr_input_result.RData")

fa_input_result = fasta_input(options$fa_file, shapeind = "all")

#save(fa_input_result, file="fa_input_result.RData")
#
##load("fa_input_result.RData")
##print(fa_input_result$gr_file)
##print(fa_input_result$dna_shape_data)
##print(fa_input_result$dna_shape_data[,,1])

dna_shape_data = fa_input_result

record_number = dim(dna_shape_data)[3]
if (record_number < 50) {
    peak_count = record_number
} else {
    peak_count = 50
}

EM_Gibbs_result = SMEM_Gibbs(
  dna_shape_data,
  # output file prefix
  filename = "EM_Gibbs",
  # number of peaks per batch
  peakCount = peak_count,
  motifLength = MOTIF_LEN,
  # number of motifs to be discovered for each peak
  motifCount = 1,
  # number of repetitions per batch
  replicates = 50,
  tolerance = 1e-06
)

#save(EM_Gibbs_result, file=paste0(prefix, "_EM_Gibbs_result.RData"))
#load("EM_Gibbs_result.RData")

# output is "A list of motif location arrays"
#print(EM_Gibbs_result)

#gr_file = fa_input_result$gr_file

motif_location_array = EM_Gibbs_result

fa_input_merge_result = fasta_motif_merge(motif_location_array,
                                         options$fa_file,
                                         motifLength = MOTIF_LEN,
                                         filename = paste0(prefix, "_peak_fa_mul"))
#print(fa_input_merge_result)
save(fa_input_merge_result, file=paste0(prefix, "_fa_input_merge_result.RData"))

#load("only_peaks_train_main_fa_input_merge_result.RData")

# motif_data_frame is shape (607, 9).
# the 9 columns are:
#   1. peakNo.
#   2. seriesNo.
#   3. chrNo.
#   4. location
#   5. genoStart
#   6. genoEnd
#   7. motiflength
#   8. mergerMotifLength
#   9. motifFrequency
#print(head(fa_input_merge_result$motif_data_frame, 8))
#print(dim(fa_input_merge_result$motif_data_frame))
#print(names(fa_input_merge_result$motif_data_frame))

# motif_shape_data is a list of length 98.
# unclear exactly what it represents. Each element of the list is, itself a list
# those sub-lists are variable in length, usually 5-7 elements though
# Each sub-list is a matrix of shape (14, x), where x is variable and must
#   relate to the length of a piece of DNA, and the 14 is for the number of shapes
#print(fa_input_merge_result$motif_shape_data[0:2])
#print(length(fa_input_merge_result$motif_shape_data))

motif_data_frame = fa_input_merge_result$motif_data_frame

motif_shape_data = fa_input_merge_result$motif_shape_data

motif_vis_pearson = motif_visualizer(
  motif_data_frame,
  motif_shape_data,
  motifFrequency = 0.5,
  shapeIndex = "all",
  motifLength = MOTIF_LEN,
  align_method = "pearson",
  ref_number = 1,
  input_form = 0,
  continue_alignment = 1,
  top_motif = 25,
  threshold = 0.6
)


#seq_list_set = motif_vis_pearson$align_motif_sequence

align_df = motif_vis_pearson$align_motif_dataframe

#shape_dfm_set = align_df$motif_shape_df_melt
#

area_dfm_set = align_df$motif_area_df_melt

save(area_dfm_set, file=paste0(prefix, "_area_dfm_set.RData"))

#
#print(seq_list_set)
#print(dna_shape_data)

#save.image()

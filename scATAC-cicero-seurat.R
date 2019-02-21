#!/usr/bin/env Rscript


#################################################################################
  #   Wrapper code adapted from:
  #     Title: Analysis of 10x Genomics Chromium single cell RNA-seq data using Seurat (version 3.0) starting with Cell Ranger output.
  #     Author: Igor Dolgalev
  #     Date: February 16th, 2019
  #     Availability: https://github.com/igordot/genomics/blob/master/scripts/scrna-10x-seurat-3.R
################################################################################

#################################################################################
#   scATAC code adapted from Cicero Docs:
#     Title: Cicero
#     Author: Hannah A. Pliner, Jay Shendure & Cole Trapnell et. al
#     Date: 2018
#     Availability: https://github.com/cole-trapnell-lab/cicero-release
#     citation("cicero")
################################################################################


#################################################################################
#   scRNAseq from Seurat:
#     Title: Seurat
#     Author: Butler et al.
#     Date: 2018
#     Availability: https://github.com/satijalab/seurat
#     citation("Seurat")
################################################################################


"
Analysis of 10x Genomics Chromium single cell ATAC-seq data using Cicero (version 1.0.14) starting with Cell Ranger output.
Basic workflow steps:
  1 - create - import counts matrix, perform initial QC, output RDS with CellDataSet,  co-acessibility, cis-co-assibility, unnorm_gene_activity, and norm_gene_activity
  2 - to_seurat - concvert to seurat object, visualize, cluster
  3 - identify - identify clusters from gene_activity based on specified clustering/resolution (higher resolution for more clusters)
Optional steps:
  normalize - gene activities of several samples
  combine - merge multiple seurat objects
  integrate - perform integration (batch correction) across multiple seurat objects
  de - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat-3.R create <analysis_dir> <sample_name> <data_dir> <gtf_dir> [--min_genes=<n> --max_genes=<n> --mt=<n>]
  scrna-10x-seurat-3.R cluster <analysis_dir> <num_dim>
  scrna-10x-seurat-3.R identify <analysis_dir> <resolution>
  scrna-10x-seurat-3.R combine <analysis_dir> <sample_analysis_dir>...
  scrna-10x-seurat-3.R integrate <analysis_dir> <num_dim> <batch_analysis_dir>...
  scrna-10x-seurat-3.R de <analysis_dir> <resolution>
  scrna-10x-seurat-3.R --help
Options:
  --min_fragments=<n> cutoff for minimum number of fragments per cell (5,000 if not specified)
  --max_fragments=<n> cutoff for maximum number of fragments per cell (98th percentile if not specified)
  --min_genes=<n>   cutoff for minimum number of genes per cell (2nd percentile if not specified)
  --max_genes=<n>   cutoff for maximum number of genes per cell (98th percentile if not specified)
  -h, --help        show this screen
" -> doc



load_libraries <- function(){

  message("\n\n ========== load libraries ========== \n\n")

  suppressPackageStartupMessages({
  library(cicero)
  library(docopt)
  library(Gviz)
  library(tidyverse)
  library(readr)
  library(glue)
  library(future)
  library(RColorBrewer)
  library(refGenome)
  library(ggplot2)
})

}

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
# Alexis Lucattini
get_this_file <- function(){
    commandArgs() %>%
       tibble::as.tibble() %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value)
}


#ata_dir -- outer dir containing filteredpeakbcmatrix
#gtf_dir -- dir to gtf file
#sample_name -- name of sample to append
#analysis_dir -- output dir

# ========== main ==========

# output width
options(width = 120)
# print warnings as they occur
options(warn = 1)
# default type for the bitmap devices such as png (should default to "cairo")
options(bitmapType = "cairo")

# retrieve the command-line arguments
suppressPackageStartupMessages(library(docopt))
opts = docopt(doc)


# dependencies
load_libraries()

path_to_this_file <- get_this_file()[1]
source(file.path(dirname(path_to_this_file),"cicero_wrap.R"))


# evaluate R expressions asynchronously when possible (such as ScaleData)
plan("multiprocess", workers = 4)
# increase the limit of the data to be shuttled between the processes from default 500MB to 50GB
options(future.globals.maxSize = 50e9)

# global settings
#colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
#colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

# analysis info
analysis_step = "unknown"
out_dir = opts$analysis_dir

# create analysis directory if starting new analysis or exit if analysis already exists
if (opts$create || opts$combine || opts$integrate) {

  if (opts$create) analysis_step = "create"
  if (opts$combine) analysis_step = "normalize"
  if (opts$combine) analysis_step = "combine"
  if (opts$integrate) analysis_step = "integrate"

  message(glue("\n\n ========== started analysis step {analysis_step} for {out_dir} ========== \n\n"))

  if (dir.exists(out_dir)) {
    stop(glue("output analysis dir {out_dir} already exists"))
  } else {
    dir.create(out_dir)
  }

  #sink(paste(opts$analysis_dir, "/", opts$sample_name, 'log.txt', sep = ""), type = "message")

  # original working dir (before moving to analysis dir)
  original_wd = getwd()

}

# set analysis directory as working directory
if (dir.exists(out_dir)) {
  setwd(out_dir)
} else {
  stop(glue("output analysis dir {out_dir} does not exist"))
}

if (opts$create) {

  # log to file
  message(glue("analysis: {out_dir}"), file = "create.log", append = TRUE)
  message(glue("seurat version: {packageVersion('cicero')}"),
    file = "create.log", append = TRUE)

  # read in 10x
  # return indata, peakinfo, and cellinfo as a list
  indata_peakinfo_cellinfo  <- read_10x(opts$data_dir)

  # get fragments per cell from indata, and output
  # return fragments_per_cell, min_cutoff and max_cutoff as a list
  fragments_per_cell <- get_fragments_per_cell(indata_peakinfo_cellinfo$indata,
    min_fragments = NULL, max_fragments = NULL)

  # plots the fragments per cell and the cutoff values
  plot_fragments_per_cell(fragments_per_cell$fragments_per_cell,
    opts$analysis_dir, opts$sample_name, fragments_per_cell$min_cutoff, fragments_per_cell$max_cutoff)

  # filter the data using the cutoff values
  # return filtered indata and filtered cellinfo
  filtered_indata_cellinfo <- filter_10x(indata_peakinfo_cellinfo$indata, fragments_per_cell$fragments_per_cell,
     indata_peakinfo_cellinfo$cellinfo, fragments_per_cell$min_cutoff,
     fragments_per_cell$max_cutoff)

  # returns celldataset
  input_cds_obj <- create_input_cds(filtered_indata_cellinfo$pass_indata,
    filtered_indata_cellinfo$pass_cellinfo, indata_peakinfo_cellinfo$peakinfo)

  # returns cicero celldataset
  cicero_cds_obj <- create_cicero_cds(input_cds_obj)

  # returns coaccessibility object
  coaccessibility <- get_coaccessibility(cicero_cds_obj)

  # returns cis coaccessibility object
  ciscoaccessibility_net <- get_ciscoaccessibility_net(coaccessibility)

  # returns unnormalized gene activities
  gene_activity <- get_gene_activity(input_cds_obj, coaccessibility, opts$gtf_dir)

  # returns normalized gene activities
  normalized_gene_activity <- norm_gene_activity(gene_activity, input_cds_obj)


  saveRDS(list(input_cds_obj, cicero_cds_obj, coaccessibility,
               ciscoaccessibility_net, gene_activity, normalized_gene_activity),
          glue("{opts$analysis_dir}/{opts$sample_name}_cicero.RDS"))
}

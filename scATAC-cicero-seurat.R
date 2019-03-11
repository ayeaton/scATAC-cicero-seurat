#!/usr/bin/env Rscript


#################################################################################
  #   Wrapper code adapted from:
  #     Title: Analysis of 10x Genomics Chromium single cell RNA-seq data using
  #            Seurat (version 3.0) starting with Cell Ranger output.
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
  2 - normalize_cicero - import unnormalized gene activities and normalize across samples
  3 - to_seurat - concvert to seurat object, visualize, cluster
  4 - identify - identify clusters from gene_activity based on specified clustering/resolution (higher resolution for more clusters)
Optional steps:
  combine - merge multiple seurat objects
  integrate - perform integration (batch correction) across multiple seurat objects
  de - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat-3.R create <analysis_dir> <sample_name> <data_dir> <gtf_dir> [--min_fragments=<n> --max_fragments=<n> --mt=<n>]
  scrna-10x-seurat-3.R normalize_cicero <analysis_dir> <analysis_name> <sample_names>... <sample_analysis_dirs>...
  scrna-10x-seurat-3.R to_seurat <analysis_dir>
  scrna-10x-seurat-3.R cicerobulk_toseurat <sample_analysis_dirs> <sample_names> <analysis_name>
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
  library(Seurat)
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

load_source_files <- function(){
  path_to_this_file <- get_this_file()[1]
  source(file.path(dirname(path_to_this_file),"cicero_wrap.R"))
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
load_source_files()

# evaluate R expressions asynchronously when possible (such as ScaleData)
plan("multiprocess", workers = 4)
# increase the limit of the data to be shuttled between the processes from default 500MB to 50GB
options(future.globals.maxSize = 50e9)


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

  file = file(paste(opts$analysis_dir, "/", opts$sample_name, 'log.txt', sep = ""),
   'w')

  sink(file, type = "message")

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
  # create RDS that contains the cicero object outputs
  cicero_obj_list <- execute_cicero_wraps(opts$analysis_dir, opts$sample_name, opts$data_dir, opts$gtf_dir, opts$min_fragments)
  saveRDS(cicero_obj_list,
          glue("{opts$analysis_dir}/{opts$sample_name}_cicero.RDS"))
} else if(opts$normalize_cicero){
  # given the location of the cicero objects, grab unnorm_gene_activity and
  # num_genes 
  list_unnormed_gene_activity <- get_bulk_gene_activity(opts$sample_names, opts$sample_analysis_dirs)
  uniform_gene_activities <- make_gene_activities_uniform(list_unnormed_gene_activity$unnorm_gene_list, opts$sample_names)
  normalize_bulk_activity <- normalize_bulk(uniform_gene_activities, 
                                            list_unnormed_gene_activity$num_genes)
  saveRDS(normalize_bulk_activity,
          glue("{opts$analysis_dir}/{opts$analysis_name}_cicero_norm.RDS"))
  save_seurat_objs(normalize_bulk_activity, out_dir,list_unnormed_gene_activity$fragment_counts)
}else if(opts$cicerobulk_toseurat){
  combine <- combine_seurat_obj(opts$sample_analysis_dirs, opts$sample_names, opt$analysis_name)
}else if(opts$to_seurat){
  cicero_obj <- read_cicero_obj(glue("{opts$analysis_dir}/{opts$sample_name}_cicero.RDS"))
  seurat_obj <- variance_plots(cicero_obj$normalized_gene_activity[[1]], opts$sample_name, out_dir)
}

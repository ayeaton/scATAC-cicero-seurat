#!/usr/bin/env Rscript

"
Analysis of 10x Genomics Chromium single cell ATAC-seq data using Cicero (version 1.0.14) starting with Cell Ranger output.
Basic workflow steps:
  1 - create - import counts matrix, perform initial QC, output RDS with CellDataSet, co-acessibility, cis-co-assibility, unnorm_gene_activity, and norm_gene_activity
  2 - to_seurat - concvert to seurat object, visualize, cluster
  3 - identify - identify clusters from gene_activity based on specified clustering/resolution (higher resolution for more clusters)
Optional steps:
  normalize - gene activities of several samples
  combine - merge multiple seurat objects
  integrate - perform integration (batch correction) across multiple seurat objects
  de - differential expression between samples/libraries within clusters

Usage:
  scrna-10x-seurat-3.R create <analysis_dir> <sample_name> <sample_dir> [--min_genes=<n> --max_genes=<n> --mt=<n>]
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


library(cicero)
library(docopt)
library(Gviz)
opts = docopt(doc)

opts$

data_dir -- outer dir containing filteredpeakbcmatrix
gtf_dir -- dir to gtf file
sample_name -- name of sample to append
analysis_dir -- output dir

# Create_CDS --------------------------------------------------------------


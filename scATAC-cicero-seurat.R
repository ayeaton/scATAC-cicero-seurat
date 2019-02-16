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

create_input_cds <- function(data_dir){
  # read in matrix data using the Matrix package
  indata <- Matrix::readMM(paste(data_dir, "filtered_peak_bc_matrix/matrix.mtx", sep = "/") )
  # binarize the matrix
  indata@x[indata@x > 0] <- 1
  
  # format cell info
  cellinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/barcodes.tsv", sep = "/") )
  row.names(cellinfo) <- cellinfo$V1
  names(cellinfo) <- "cells"
  
  # format peak info
  peakinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/peaks.bed", sep = "/") )
  names(peakinfo) <- c("chr", "bp1", "bp2")
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name
  
  row.names(indata) <- row.names(peakinfo)
  colnames(indata) <- row.names(cellinfo)
  
  # make CDS
  fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
  pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
  input_cds <-  suppressWarnings(newCellDataSet(indata,
                                                phenoData = pd,
                                                featureData = fd,
                                                expressionFamily=VGAM::binomialff(),
                                                lowerDetectionLimit=0))
  input_cds@expressionFamily@vfamily <- "binomialff"
  input_cds <- monocle::detectGenes(input_cds)
  
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
  
  set.seed(2017)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  # *** if you are using Monocle 3, you need to run the following line as well!
  #input_cds <- preprocessCDS(input_cds, norm_method = "none")
  input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                               reduction_method = 'tSNE', norm_method = "none")
  return(input_cds)
}

create_cicero_cds <- function(input_cds){
  tsne_coords <- t(reducedDimA(input_cds))
  row.names(tsne_coords) <- row.names(pData(input_cds))
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
  return(cicero_cds)
}

read_gtf <- function(gtf_dir){  
  ens <- ensemblGenome()
  read.gtf <- read.gtf(ens, gft_dir))
  gene_table <- getGeneTable(ens)
  return(gene_table)
}

create_coaccessibility_obj <- function(cicero_cds, input_cds){

  cleaned_gene_table <- data.frame(chromosome = paste("chr",gene_table$seqid, sep = ""),
                                  start = gene_table$start, 
                                  end = gene_table$end,
                                  gene = gene_table$gene_name)
  
  cleaned_gene_table <- cleaned_gene_table[!is.na(cleaned_gene_table$transcript),]
  
  coaccessibility <- run_cicero(cicero_cds, cleaned_gene_table) # Takes a few minutes to run
}


  cis_coaccessibility_network <- generate_ccans(coaccessibility)
}




input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
unnorm_gene_activity <- build_gene_activity_matrix(input_cds, coaccessibility)

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
norm_gene_activities <- normalize_gene_activities(unnorm_gene_activity, num_genes)





saveRDS(coaccessibility, cis_coaccessibility_network, unnorm_gene_activity, norm_gene_activities, paste(analysis_dir, "/", sample_name, "_cicero.RDS"))



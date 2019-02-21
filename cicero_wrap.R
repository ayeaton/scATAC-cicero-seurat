#!/usr/bin/env Rscript
read_10x <- function(data_dir){
  # read in matrix data using the Matrix package
  indata <- Matrix::readMM(paste(data_dir, "filtered_peak_bc_matrix/matrix.mtx", sep = "/") )
  # binarize the matrix
  indata@x[indata@x > 0] <- 1

  # format cell info
  cellinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/barcodes.tsv", sep = "/"), stringsAsFactors = F)
  row.names(cellinfo) <- cellinfo$V1
  names(cellinfo) <- "cells"

  # format peak info
  peakinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/peaks.bed", sep = "/") )
  names(peakinfo) <- c("chr", "bp1", "bp2")
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name

  row.names(indata) <- row.names(peakinfo)
  colnames(indata) <- row.names(cellinfo)
  return(list(indata=indata, peakinfo=peakinfo, cellinfo= cellinfo))
}

get_fragments_per_cell <- function(indata, min_fragments = NULL, max_fragments = NULL){
  fragments_per_cell <- data.frame(fragment_count=Matrix::colSums(indata))
  # get the quantile vals
  low_quantile <- fragments_per_cell$fragment_count %>% quantile(0.02) %>% round(1)
  high_quantile <- fragments_per_cell$fragment_count %>% quantile(0.98) %>% round(1)

  # set cutoff values
  if(is.null(min_fragments)){
    min_cutoff = low_quantile
  } else{
    min_cutoff = min_fragments
    names(min_cutoff) <- "min_fragments"
  }

  if(is.null(max_fragments)){
    max_cutoff = high_quantile
  } else{
    max_cutoff = max_fragments
    names(max_cutoff) <- "max_framents"
  }

  return(list(fragments_per_cell=fragments_per_cell, min_cutoff=min_cutoff, max_cutoff=max_cutoff))
}

plot_fragments_per_cell <- function(fragments_per_cell, analysis_dir, sample_name, min_cutoff, max_cutoff){

  message(names(fragments_per_cell) )
  boxplot_colors <- c("orange", "sky blue")


  fragments_per_cell$cutoff_col <- rep("cut", nrow(fragments_per_cell))
  pass_ind <- which((fragments_per_cell > min_cutoff) & (fragments_per_cell < max_cutoff))
  fragments_per_cell$cutoff_col[pass_ind] <- "kept"

  boxplot_raw <- ggplot(fragments_per_cell, aes(x = "fragment_count" , y=fragment_count)) + geom_boxplot() +
    geom_jitter(aes(colour = cutoff_col), shape=16,
   position=position_jitter(0.2)) + theme_minimal() +
   scale_color_manual(values = boxplot_colors)

  ggsave(boxplot_raw,    filename=glue("{analysis_dir}/{sample_name}_fragments_per_cell.png"), device = "png")



  fragments_per_cell$log10_fragment_count <- log10(fragments_per_cell$fragment_count)

  log10_boxplot <- ggplot(fragments_per_cell, aes(x = "fragment_count" , y=log10_fragment_count)) + geom_boxplot() +
    geom_jitter(aes(colour = cutoff_col), shape=16,
   position=position_jitter(0.2)) + theme_minimal() +
   scale_color_manual(values = boxplot_colors)

  ggsave(log10_boxplot, filename=glue("{analysis_dir}/{sample_name}_log10_fragments_per_cell.png"), device = "png")
}

# FiLTER Cell info
filter_10x <- function(indata, fragments_per_cell, cellinfo, min_cutoff=NULL, max_cutoff=NULL){
  # filter
  pass_ind <- which((fragments_per_cell > min_cutoff) & (fragments_per_cell < max_cutoff))

  pass_cells <- indata[,pass_ind]

  pass_cellinfo <- as.data.frame(cellinfo[pass_ind,])

  names(pass_cellinfo) <- "cells"
  rownames(pass_cellinfo) <- pass_cellinfo$cells
  return(list(pass_indata=pass_cells, pass_cellinfo=pass_cellinfo))

}

create_input_cds <- function(indata, cellinfo, peakinfo){
  #read in 10x data
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

get_gene_table <- function(gtf_dir){
  setwd(dirname(gtf_dir))
  ens <- ensemblGenome()
  read.gtf <- read.gtf(ens, basename(gtf_dir))
  gene_table <- getGeneTable(ens)
  cleaned_gene_table <- data.frame(chromosome = paste("chr",gene_table$seqid, sep = ""),
                                   start = gene_table$start,
                                   end = gene_table$end,
                                   gene = gene_table$gene_name)
  cleaned_gene_table <- cleaned_gene_table[!is.na(cleaned_gene_table$gene),]
  return(cleaned_gene_table)
}

#ensemble is taking forever to respond. Re do this with the genome file from Ensemble
get_chrom_size <- function(){
  mouse.mm10.genome <- data.frame(V1 = c(paste("chr", 1:19, sep = ""), "chrY", "chrX"),
                                  V2 <- c(195471971, 182113224, 160039680, 156508116,
                                          151834684, 149736546, 145441459, 129401213,
                                          124595110, 130694993, 122082543, 120129022,
                                          120421639, 	124902244, 104043685, 98207768,
                                          94987271, 90702639, 61431566, 171031299,
                                          91744698))
}

get_coaccessibility <- function(cicero_cds){
  chrom_size <- get_chrom_size()
  coaccessibility <- run_cicero(cicero_cds, chrom_size) # Takes a few minutes to run
  return(coaccessibility)
}

get_ciscoaccessibility_net <- function(coaccessibility){
  cis_coaccessibility_network <- generate_ccans(coaccessibility)
  return(cis_coaccessibility_network)
}

get_gene_activity <- function(input_cds, coaccessibility, gtf_dir, verbose = FALSE){
  say <- message
  if (!verbose){
    say <- function(message) invisible(NULL)
  }
  wd <- getwd()
  gene_table <- get_gene_table(gtf_dir)
  setwd(wd)
  input_cds <- annotate_cds_by_site(input_cds, gene_table)
  unnorm_gene_activity <- build_gene_activity_matrix(input_cds, coaccessibility)
}

norm_gene_activity <- function(unnorm_gene_activity, input_cds){
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  # normalize
  norm_gene_activities <- normalize_gene_activities(unnorm_gene_activity, num_genes)
}

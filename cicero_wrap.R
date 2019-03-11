#!/usr/bin/env Rscript
read_10x <- function(data_dir){
  # read in matrix data using the Matrix package
  indata <- Matrix::readMM(paste(data_dir, "filtered_peak_bc_matrix/matrix.mtx",
    sep = "/") )
  # binarize the matrix
  indata@x[indata@x > 0] <- 1

  # format cell info
  cellinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/barcodes.tsv",
    sep = "/"), stringsAsFactors = F)
  row.names(cellinfo) <- cellinfo$V1
  names(cellinfo) <- "cells"

  # format peak info
  peakinfo <- read.table(paste(data_dir,"filtered_peak_bc_matrix/peaks.bed",
    sep = "/") )
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
    min_cutoff = as.numeric(min_fragments)
    names(min_cutoff) <- "min_fragments"
  }

  if(is.null(max_fragments)){
    max_cutoff = high_quantile
  } else{
    max_cutoff = as.numeric(max_fragments)
    names(max_cutoff) <- "max_framents"
  }

  return(list(fragments_per_cell=fragments_per_cell, min_cutoff=min_cutoff,
    max_cutoff=max_cutoff))
}

plot_fragments_per_cell <- function(fragments_per_cell, analysis_dir,
  sample_name, min_cutoff, max_cutoff){

  message(names(fragments_per_cell) )
  boxplot_colors <- c("orange", "sky blue")

  print(min_cutoff)
  print(max_cutoff)
  fragments_per_cell$cutoff_col <- rep("cut", nrow(fragments_per_cell))
  pass_ind <- which((fragments_per_cell > min_cutoff) & (fragments_per_cell < max_cutoff))
  fragments_per_cell$cutoff_col[pass_ind] <- "kept"

  boxplot_raw <- ggplot(fragments_per_cell, aes(x = "fragment_count" , y=fragment_count)) +
    geom_boxplot() + geom_jitter(aes(colour = cutoff_col), shape=16, position=position_jitter(0.2)) +
    theme_minimal() + scale_color_manual(values = boxplot_colors)

  ggsave(boxplot_raw,    filename=glue("{analysis_dir}/{sample_name}_fragments_per_cell.png"),
    device = "png")



  fragments_per_cell$log10_fragment_count <- log10(fragments_per_cell$fragment_count)

  log10_boxplot <- ggplot(fragments_per_cell, aes(x = "fragment_count" , y=log10_fragment_count)) +
    geom_boxplot() + geom_jitter(aes(colour = cutoff_col), shape=16, position=position_jitter(0.2)) +
    theme_minimal() + scale_color_manual(values = boxplot_colors)

  ggsave(log10_boxplot, filename=glue("{analysis_dir}/{sample_name}_log10_fragments_per_cell.png"),
    device = "png")
}

# FiLTER Cell info
filter_10x <- function(indata, fragments_per_cell, cellinfo, min_cutoff=NULL, max_cutoff=NULL){
  # filter
  pass_ind <- which((fragments_per_cell > min_cutoff) & (fragments_per_cell < max_cutoff))

  pass_cells <- indata[,pass_ind]

  pass_cellinfo <- as.data.frame(cellinfo[pass_ind,])
  
  pass_fragments_per_cell <- fragments_per_cell[pass_ind,]

  names(pass_cellinfo) <- "cells"
  rownames(pass_cellinfo) <- pass_cellinfo$cells
  return(list(pass_indata=pass_cells, pass_cellinfo=pass_cellinfo, fragments_per_cell = pass_fragments_per_cell))

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
  saveRDS(file = "cleaned_gene_table.RDS", cleaned_gene_table)
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
  return(list(norm_gene_activities = norm_gene_activities, num_genes = num_genes))
}

get_bulk_gene_activity <- function(sample_names, sample_analysis_dirs){

  if (length(sample_names) != length(sample_analysis_dirs)){
    stop(glue("number of sample names and sample analysis directories is not equal"))
  }
  
  unnorm_gene_list = list()
  num_genes = list()
  fragment_counts = list()
  for (i in 1:length(sample_analysis_dirs)) {
    sample_name <- sample_names[i]
    sample_analysis_dir = sample_analysis_dirs[i]
    sample_cicero_rds = glue("{sample_analysis_dir}/{sample_name}_cicero.RDS")

    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir[[1]])) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_cicero_rds[[1]])) stop(glue("seurat object rds {sample_cicero_rds} does not exist"))

    # load cicero object
     tmp = readRDS(sample_cicero_rds)
     unnorm_gene_list[[i]] = tmp$gene_activity
     num_genes[[i]] = tmp$num_genes
     fragment_counts[[i]] = tmp$fragments_per_cell
     rm(tmp)
  }
  return(list(unnorm_gene_list = unnorm_gene_list, num_genes = num_genes, fragment_counts = fragment_counts))
}

make_gene_activities_uniform <- function(list_gene_activity, sample_names){
  all_genes <- lapply(list_gene_activity, function(x) rownames(x))
  all_genes <- unique(unlist(all_genes))
  all_genes <- as.data.frame(all_genes, stringsAsFactors = F)
  rownames(all_genes) <- all_genes[,1]
  
  uniform_gene_activity <- lapply(list_gene_activity, function(x){
    tmp <- merge(all_genes, as.data.frame(x), by = "row.names", all = T)
    NA_to_zero <- tmp %>% 
      mutate_all(funs(replace(., is.na(.), 0)))
    rownames(NA_to_zero) <- NA_to_zero$Row.names
    NA_to_zero <- NA_to_zero[,-c(1,2)]
    out <- Matrix(as.matrix(NA_to_zero), sparse = TRUE)
  })
  names(uniform_gene_activity) <- sample_names
  return(uniform_gene_activity)
}

normalize_bulk <- function(list_gene_activity, num_genes){
  num_genes <- unlist(num_genes)
  bulk_norm <- normalize_gene_activities(list_gene_activity, num_genes)
  return(bulk_norm)
}

save_seurat_objs <- function(normalize_bulk_activity, out_dir, fragment_counts){
  for (i in 1:length(normalize_bulk_activity)){
    sub_dir = glue("{out_dir}/{names(normalize_bulk_activity)[i]}")
    print(sub_dir)
    #make new directory for each s obj
     if (dir.exists(sub_dir)) {
      stop(glue("output analysis dir {sub_dir} already exists"))
     } else {
      dir.create(sub_dir)
     }
    out <- variance_plots(normalize_bulk_activity[[i]], names(normalize_bulk_activity)[i], sub_dir, fragment_counts[[i]], TRUE)
  }
}

execute_cicero_wraps <- function(analysis_dir, sample_name, data_dir, gtf_dir, min_fragments){

  # read in 10x
  # return indata, peakinfo, and cellinfo as a list
  indata_peakinfo_cellinfo  <- read_10x(data_dir)

  # get fragments per cell from indata, and output
  # return fragments_per_cell, min_cutoff and max_cutoff as a list
  fragments_per_cell <- get_fragments_per_cell(indata_peakinfo_cellinfo$indata,
    min_fragments = min_fragments, max_fragments = NULL)

  # plots the fragments per cell and the cutoff values
  plot_fragments_per_cell(fragments_per_cell$fragments_per_cell,
    analysis_dir, sample_name, fragments_per_cell$min_cutoff, fragments_per_cell$max_cutoff)

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
  gene_activity <- get_gene_activity(input_cds_obj, coaccessibility, gtf_dir)

  # returns normalized gene activities
  normalized_gene_activity <- norm_gene_activity(gene_activity, input_cds_obj)

  return(list(input_cds_obj = input_cds_obj, cicero_cds_obj = cicero_cds_obj,
    coaccessibility = coaccessibility,ciscoaccessibility_net = ciscoaccessibility_net,
    gene_activity = gene_activity,
    normalized_gene_activity = normalized_gene_activity$norm_gene_activities,
    num_genes = normalized_gene_activity$num_genes,
    fragments_per_cell = filtered_indata_cellinfo$fragments_per_cell))
}

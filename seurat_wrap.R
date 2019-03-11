#!/usr/bin/env Rscript

read_cicero_obj <- function(cicero_obj){
  cicero_rds <- readRDS(cicero_obj)
}
variance_plots <- function(gene_activity, sample_name, out_dir, fragment_counts, plot = FALSE){
  message("\n\n ========== Seurat::FindVariableFeatures()========== \n\n")

  print(length(fragment_counts))
  
  seurat_obj <- CreateSeuratObject(gene_activity, project = sample_name)
  seurat_obj$counts <- fragment_counts
  
  seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst",
    verbose = FALSE)

  if(plot){
    var_plot = VariableFeaturePlot(seurat_obj, pt.size = 0.5)
    ggsave(glue("{out_dir}/variance.features.png"), plot = var_plot, width = 12, height = 5,
      units = "in")
  }

  message("\n\n ========== Seurat::ScaleData()========== \n\n")

  seurat_obj = ScaleData(seurat_obj, verbose = FALSE)

  message("\n\n ========== Seurat::RunPCA()========== \n\n")
  # PCA on the scaled data
  # PCA calculation stored in object[["pca"]]
  num_pcs = 50
  if (ncol(seurat_obj) < 100) num_pcs = 20
  if (ncol(seurat_obj) < 25) num_pcs = 5
  seurat_obj = RunPCA(seurat_obj, assay = "RNA", features = VariableFeatures(seurat_obj),
    npcs = num_pcs, verbose = FALSE)

  if(plot){
    # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
    pca_plot =
      DimPlot(
        seurat_obj, cells = sample(colnames(seurat_obj)), group.by = "orig.ident", reduction = "pca",
        pt.size = 0.5
      ) +
      theme(aspect.ratio = 1)
    ggsave(glue("{out_dir}/variance.pca.png"), plot = pca_plot, width = 8, height = 6, units = "in")
  }

  message("\n\n ========== Seurat::DimHeatmap() ========== \n\n")

  if(plot){
    # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
    if (num_pcs > 15) {
      png(glue("{out_dir}/variance.pca.heatmap.png"), res = 300, width = 10, height = 16, units = "in")
        DimHeatmap(seurat_obj, reduction = "pca", dims = 1:15, nfeatures = 20, cells = 250, fast = TRUE)
      dev.off()
    }
  }

  message("\n\n ========== Seurat::PCElbowPlot() ========== \n\n")

  if(plot){
    # a more ad hoc method for determining PCs to use, draw cutoff where there is a clear elbow in the graph
    elbow_plot = ElbowPlot(seurat_obj, reduction = "pca", ndims = num_pcs)
    ggsave(glue("{out_dir}/variance.pca.elbow.png"), plot = elbow_plot, width = 8, height = 5, units = "in")

    # resampling test inspired by the jackStraw procedure - very slow, so skip for large projects (>10,000 cells)
    if (ncol(seurat_obj) < 10000) {

      message("\n\n ========== Seurat::JackStraw() ========== \n\n")

      # determine statistical significance of PCA scores
      seurat_obj = JackStraw(seurat_obj, assay = "RNA", reduction = "pca", dims = num_pcs, verbose = FALSE)

      # compute Jackstraw scores significance
      seurat_obj = ScoreJackStraw(seurat_obj, reduction = "pca", dims = 1:num_pcs, do.plot = FALSE)

      # plot the results of the JackStraw analysis for PCA significance
      # significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line)
      jackstraw_plot =
        JackStrawPlot(seurat_obj, reduction = "pca", dims = 1:num_pcs) +
        guides(col = guide_legend(ncol = 2))
      ggsave(glue("{out_dir}/variance.pca.jackstraw.png"), plot = jackstraw_plot, width = 12, height = 6, units = "in")
      }
  }
  saveRDS(seurat_obj, glue("{out_dir}/{sample_name}_seurat_obj.RDS"))
}

# merge multiple Seurat objects
combine_seurat_obj = function(sample_analysis_dirs, sample_names, out_dir, plot = TRUE) {
  
  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
  
  message("\n\n ========== combine samples ========== \n\n")
  
  cicero_obj_list = list()
  for (i in 1:length(sample_analysis_dirs)) {
    sample_name = sample_names[[i]]
    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{sample_analysis_dir}")
    sample_cicero_rds = glue("{sample_analysis_dir}/{sample_name}_seurat_obj.RDS")
    
    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_cicero_rds)) stop(glue("seurat object rds {sample_cicero_rds} does not exist"))
    
    # load seurat object
    cicero_obj_list[[i]] = readRDS(sample_cicero_rds)
    
    # clean up object
    cicero_obj_list[[i]]@assays$RNA@var.features = vector()
    cicero_obj_list[[i]]@assays$RNA@scale.data = matrix()
    cicero_obj_list[[i]]@reductions = list()
    cicero_obj_list[[i]]@meta.data = cicero_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))
    
    # print single sample sample stats
    # sample_name = seurat_obj_list[[i]]@meta.data[1, "orig.ident"] %>% as.character()
    sample_name = cicero_obj_list[[i]]$orig.ident[1] %>% as.character()

  }
  
  # merge
  seurat_obj = merge(cicero_obj_list[[1]], cicero_obj_list[2:length(cicero_obj_list)])
  rm(cicero_obj_list)
  

  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(seurat_obj, assay = "RNA", slot = "counts") > 0)
  min_cells = 10
  if (ncol(seurat_obj) > 100000) { min_cells = 50 }
  filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
  seurat_obj = subset(seurat_obj, features = filtered_genes)
  

  # print gene/cell minimum cutoffs
  min_cells = Matrix::rowSums(GetAssayData(seurat_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  min_genes = Matrix::colSums(GetAssayData(seurat_obj, assay = "RNA", slot = "counts") > 0) %>% min()


  # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
  num_matrix_elements = GetAssayData(seurat_obj, assay = "RNA", slot = "counts") %>% length()
  if (num_matrix_elements < 2^31) {
    
    # save raw counts matrix
    counts_raw = GetAssayData(seurat_obj, assay = "RNA", slot = "counts") %>% as.matrix()
    counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    # write_csv(counts_raw, path = "counts.raw.csv.gz")
    # fwrite(counts_raw, file = "counts.raw.csv", sep = ",")
    # R.utils::gzip("counts.raw.csv")
    # rm(counts_raw)
    
    # save counts matrix as a basic gzipped text file
    counts_norm = GetAssayData(seurat_obj, assay = "RNA") %>% as.matrix() %>% round(3)
    counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    # write_csv(counts_norm, path = "counts.normalized.csv.gz")
    # fwrite(counts_norm, file = "counts.normalized.csv", sep = ",")
    # R.utils::gzip("counts.normalized.csv")
    # rm(counts_norm)
    
  }
  
  Sys.sleep(1)
  
  seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst",
                                    verbose = FALSE)
  
  if(plot){
    var_plot = VariableFeaturePlot(seurat_obj, pt.size = 0.5)
    ggsave(glue("{out_dir}/variance.features.png"), plot = var_plot, width = 12, height = 5,
           units = "in")
  }
  
  message("\n\n ========== Seurat::ScaleData()========== \n\n")
  
  seurat_obj = ScaleData(seurat_obj, verbose = FALSE)
  
  message("\n\n ========== Seurat::RunPCA()========== \n\n")
  # PCA on the scaled data
  # PCA calculation stored in object[["pca"]]
  num_pcs = 50
  if (ncol(seurat_obj) < 100) num_pcs = 20
  if (ncol(seurat_obj) < 25) num_pcs = 5
  seurat_obj = RunPCA(seurat_obj, assay = "RNA", features = VariableFeatures(seurat_obj),
                      npcs = num_pcs, verbose = FALSE)
  
  if(plot){
    # plot the output of PCA analysis (shuffle cells so any one group does not appear overrepresented due to ordering)
    pca_plot =
      DimPlot(
        seurat_obj, cells = sample(colnames(seurat_obj)), group.by = "orig.ident", reduction = "pca",
        pt.size = 0.5
      ) +
      theme(aspect.ratio = 1)
    ggsave(glue("{out_dir}/variance.pca.png"), plot = pca_plot, width = 8, height = 6, units = "in")
  }
  
  message("\n\n ========== Seurat::DimHeatmap() ========== \n\n")
  
  if(plot){
    # PCHeatmap (former) allows for easy exploration of the primary sources of heterogeneity in a dataset
    if (num_pcs > 15) {
      png(glue("{out_dir}/variance.pca.heatmap.png"), res = 300, width = 10, height = 16, units = "in")
      DimHeatmap(seurat_obj, reduction = "pca", dims = 1:15, nfeatures = 20, cells = 250, fast = TRUE)
      dev.off()
    }
  }
  
  saveRDS(seurat_obj, glue("{out_dir}/merged_seurat_obj.RDS"))
  
}

# integrate multiple Seurat objects
integrate_seurat_obj = function(original_wd, sample_analysis_dirs, sample_names, num_dim) {
  
  # check if the inputs seems reasonable
  if (length(sample_analysis_dirs) < 2) stop("must have at least 2 samples to merge")
  num_dim = as.integer(num_dim)
  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)
  
  message("\n\n ========== integrate samples ========== \n\n")
  
  seurat_obj_list = list()
  var_genes_list = list()
  exp_genes = c()
  for (i in 1:length(sample_analysis_dirs)) {
    
    sample_analysis_dir = sample_analysis_dirs[i]
    sample_analysis_dir = glue("{original_wd}/{sample_analysis_dir}")
    sample_seurat_rds = glue("{sample_analysis_dir}/seurat_obj.rds")
    
    # check if analysis dir is valid
    if (!dir.exists(sample_analysis_dir)) stop(glue("dir {sample_analysis_dir} does not exist"))
    # check if seurat object exists
    if (!file.exists(sample_seurat_rds)) stop(glue("seurat object rds {sample_seurat_rds} does not exist"))
    
    # load seurat object
    seurat_obj_list[[i]] = readRDS(sample_seurat_rds)
    # sample_name = seurat_obj_list[[i]]@meta.data[1, "orig.ident"] %>% as.character()
    sample_name = seurat_obj_list[[i]]$orig.ident[1] %>% as.character()
    
    # clean up object
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))
    
    # save expressed genes keeping only genes present in all the datasets (for genes to integrate in IntegrateData)
    if (length(exp_genes) > 0) {
      exp_genes = intersect(exp_genes, rownames(seurat_obj_list[[i]])) %>% sort()
    } else {
      exp_genes = rownames(seurat_obj_list[[i]])
    }
    
    # save variable genes
    var_genes_list[[sample_name]] = VariableFeatures(seurat_obj_list[[i]])
    
    # print single sample sample stats
    message(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"))
    write(glue("sample {sample_name} dir: {basename(sample_analysis_dir)}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} cells: {ncol(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"))
    write(glue("sample {sample_name} genes: {nrow(seurat_obj_list[[i]])}"), file = "create.log", append = TRUE)
    message(" ")
    
  }
  
  # euler plot of variable gene overlaps (becomes unreadable and can take days for many overlaps)
  if (length(var_genes_list) < 8) {
    colors_euler = colors_samples[1:length(var_genes_list)]
    euler_fit = euler(var_genes_list, shape = "ellipse")
    euler_plot = plot(euler_fit,
                      fills = list(fill = colors_euler, alpha = 0.7),
                      edges = list(col = colors_euler))
    png("variance.vargenes.euler.png", res = 200, width = 5, height = 5, units = "in")
    print(euler_plot)
    dev.off()
  }
  
  # upset plot of variable gene overlaps
  png("variance.vargenes.upset.png", res = 200, width = 8, height = 5, units = "in")
  upset(fromList(var_genes_list), nsets = 50, nintersects = 15, order.by = "freq", mb.ratio = c(0.5, 0.5))
  dev.off()
  
  message("\n\n ========== Seurat::FindIntegrationAnchors() ========== \n\n")
  
  # find the integration anchors
  anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000, dims = 1:num_dim)
  rm(seurat_obj_list)
  
  message("\n\n ========== Seurat::IntegrateData() ========== \n\n")
  
  # integrating all genes may cause issues and may not add any relevant information
  # integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim, features.to.integrate = exp_genes)
  integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim)
  rm(anchors)
  
  # after running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix
  # the original (uncorrected values) are still stored in the object in the “RNA” assay
  
  # switch to integrated assay
  DefaultAssay(integrated_obj) = "integrated"
  
  # print integrated sample stats
  message(glue("integrated unfiltered cells: {ncol(integrated_obj)}"))
  write(glue("integrated unfiltered cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated unfiltered genes: {nrow(integrated_obj)}"))
  write(glue("integrated unfiltered genes: {nrow(integrated_obj)}"), file = "create.log", append = TRUE)
  
  # filter poorly expressed genes (detected in less than 10 cells)
  filtered_genes = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0)
  min_cells = 10
  if (ncol(integrated_obj) > 100000) { min_cells = 50 }
  filtered_genes = filtered_genes[filtered_genes >= min_cells] %>% names() %>% sort()
  integrated_obj = subset(integrated_obj, features = filtered_genes)
  
  # print integrated sample stats
  message(glue("integrated cells: {ncol(integrated_obj)}"))
  write(glue("integrated cells: {ncol(integrated_obj)}"), file = "create.log", append = TRUE)
  message(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"))
  write(glue("integrated genes: {nrow(GetAssayData(integrated_obj, assay = 'RNA'))}"), file = "create.log", append = TRUE)
  
  # print gene/cell minumum cutoffs
  min_cells = Matrix::rowSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  min_genes = Matrix::colSums(GetAssayData(integrated_obj, assay = "RNA", slot = "counts") > 0) %>% min()
  message(glue("min cells per gene: {min_cells}"))
  write(glue("min cells per gene: {min_cells}"), file = "create.log", append = TRUE)
  message(glue("min genes per cell: {min_genes}"))
  write(glue("min genes per cell: {min_genes}"), file = "create.log", append = TRUE)
  
  # check that the full counts table is small enough to fit into an R matrix (max around 100k x 21k)
  num_matrix_elements = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% length()
  if (num_matrix_elements < 2^31) {
    
    # save raw counts matrix
    counts_raw = GetAssayData(integrated_obj, assay = "RNA", slot = "counts") %>% as.matrix()
    counts_raw = counts_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    # write_csv(counts_raw, path = "counts.raw.csv.gz")
    fwrite(counts_raw, file = "counts.raw.csv", sep = ",")
    R.utils::gzip("counts.raw.csv")
    rm(counts_raw)
    
    # save normalized counts matrix
    counts_norm = GetAssayData(integrated_obj, assay = "RNA") %>% as.matrix() %>% round(3)
    counts_norm = counts_norm %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    # write_csv(counts_norm, path = "counts.normalized.csv.gz")
    fwrite(counts_norm, file = "counts.normalized.csv", sep = ",")
    R.utils::gzip("counts.normalized.csv")
    rm(counts_norm)
    
    # save integrated counts matrix
    # counts_int = GetAssayData(integrated_obj, assay = "integrated") %>% as.matrix() %>% round(3)
    # counts_int = counts_int %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
    # write_csv(counts_int, path = "counts.integrated.csv.gz")
    # rm(counts_int)
    
  }
}
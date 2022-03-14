#' basicSeurat
#' 
#' Process one sample with seurat. Blindly use the provided parameters (for quick visualization only!)
#' Useful to create diagnosis graphs but not for fine exploration of the results.
#' @param ssce, Seurat Object
#' @param npcs, number of PCA values to use
#' @param nfeatures, number of features to use
#' @param SCT, logical, use SCTransform
#' @param var.to.regress, values to regress
#' @param resolution, FindCluster Resolution
#' @param rv.th, if nfeature == NULL, used residual variance cutoff instead
#' @param scale.factor, if SCT=FALSE which value to scale norm
basicSeurat <- function(ssce = NULL, nfeatures = 3000, rv.th = NULL, npcs = 25, SCT = TRUE, scale.factor = 10000, resolution = 0.5, var.to.regress = FALSE){
  
  if(SCT == FALSE){
    ssce <- NormalizeData(ssce, scale.factor = scale.factor)
    ssce <- FindVariableFeatures(ssce, selection.method = "vst", nfeatures = nfeatures)
    if(!is.null(var.to.regress )){
      ssce <- ScaleData(ssce, var.to.regress = var.to.regress)
    } else {
      ssce <- ScaleData(ssce)
    }
    ssce <- RunPCA(ssce, features = VariableFeatures(object = ssce), npcs = npcs)
    # ElbowPlot(ssce, ndims = 50)
    ssce <- FindNeighbors(ssce, dims = 1:npcs)
    ssce <- FindClusters(ssce, resolution = resolution)
    ssce <- RunUMAP(ssce, dims = 1:npcs)
  } else {
    if(!is.null(var.to.regress)){
      if(is.null(nfeatures)){
        ssce <- SCTransform(ssce, variable.features.rv.th = rv.th, variable.features.n = NULL, vars.to.regress = var.to.regress, seed.use = 42)
      } else {
        ssce <- SCTransform(ssce, vars.to.regress = var.to.regress, variable.features.n = nfeatures, seed.use = 42)
      }  
    } else {
      if(is.null(nfeatures)){
        ssce <- SCTransform(ssce, variable.features.rv.th = rv.th, variable.features.n = NULL, seed.use = 42)
      } else {
        ssce <- SCTransform(ssce, variable.features.n = nfeatures, seed.use = 42)
      }
    }
    ssce <- RunPCA(ssce)
    ssce <- FindNeighbors(ssce, dims = 1:npcs)
    ssce <- FindClusters(ssce, resolution = resolution)
    ssce <- RunUMAP(ssce, dims = 1:npcs, n.neighbors = 10, min.dist = 0.1, umap.method = "uwot", return.model=TRUE)
  }
  
  return(ssce)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

## ARCHIVE
#### Not used in the manuscripts:

#' read10xData
#' 
#' Read 10x Data as a Seurat Object and perform the first filters
#' Can flag mitochondrial genes, remove specific genes, add a GROUP column and flag droplets.
#' @param gex.path, cellRanger "filtered_feature_bc_matrix" folder.
#' @param stat.path, cellRanger "outs/metrics_summary.csv" file.
#' @param doublets.path, "*_filtered_cellIDs.txt", cell names NOT considered as doublets. if present, will flag doublets.
#' @param sample.id, ID of the sample (will be concatenated to the cell_ids)
#' @param mitochondrial.gene, vector of mitochondrial gene IDs
#' @param excluded.genes, vecor of gene IDs that should be filtered out (e.g., MALAT1, RPL*, ...)
#' @param GROUP.function, case_when wrapper function to provide the "GROUP" metadata info
read10xData <- function(gex.path = NULL, stat.path = NULL, doublets.path = NULL, sample.id = NULL, mitochondrial.gene = NULL, excluded.genes = NULL, GROUP.function = NULL){
  
  message(paste0("Reading: ", sample.id))
  
  # 2.1. Read the 10x Object & clean-up names
  expr.mtrx <- Read10X(gex.path)
  colnames(expr.mtrx) <- paste(sample.id, colnames(expr.mtrx), sep = "_")
  colnames(expr.mtrx) <- str_replace_all(colnames(expr.mtrx), pattern = "-1", replacement = "")
  sce <- CreateSeuratObject(expr.mtrx)
  
  # 2.2. Identify doublets
  if(!is.null(doublets.path)){
    single.id <- names(suppressWarnings(read.csv(doublets.path, as.is = TRUE, check.names = FALSE, quote = "'")))
    single.id <- str_replace_all(paste0(sample.id, "_", single.id), pattern = "-1", replacement = "")
  } 
  
  # 2.3. Add metadata: Mitochondrial genes
  sce <- PercentageFeatureSet(sce, features = mitochondrial.gene, col.name = "percent.mito")
  
  # 2.4. Remove some unwanted genes
  if(!is.null(excluded.genes)){
    sce <- sce[!row.names(sce) %in% excluded.genes,]
  }
  
  # 2.5. Get the sample statistics 
  statistics <- tibble(ID = sample.id,
                       orig.ident = colnames(sce),
                       nCells.total = ncol(sce),
                       nFeatures = sce$nFeature_RNA,
                       nUMI = sce$nCount_RNA,
                       mito.perc = sce$percent.mito)
  
  if(!is.null(doublets.path)){
    statistics <- mutate(statistics,
                         nCells.nonDoublet = length(single.id),
                         isDoublet = ! colnames(sce) %in% single.id,
                         singlet.perc = round(100*(nCells.nonDoublet/nCells.total),2))
  } else {
    statistics <- mutate(statistics,
                         nCells.nonDoublet = NA,
                         isDoublet = NA,
                         singlet.perc = NA)
  }
  
  if(!is.null(GROUP.function)){
    statistics <- mutate(statistics,
                         GROUP = GROUP.function(statistics$ID))
  } else {
    statistics <- mutate(statistics,
                         GROUP = sample.id)
  }
  
  # 2.6 Add the mapping statistics
  statistics <- as.data.frame(bind_cols(statistics, read_delim(stat.path, delim=",")))
  colnames(statistics) <- c("ID","orig.ident",
                            "nCells.total",
                            "nFeatures",
                            "nUMI",
                            "mito.perc",
                            "nCells.nonDoublet",
                            "isDoublet",
                            "singlet.perc",
                            "GROUP",
                            "Estimated_Number_of_Cells",
                            "Mean_Reads_per_Cell",
                            "Median_Genes_per_Cell",
                            "Number_of_Reads",
                            "Valid_Barcodes",
                            "Sequencing_Saturation",
                            "Q30_Bases_in_Barcode",
                            "Q30_Bases_in_RNA_Read",
                            "Q30_Bases_in_UMI",
                            "Reads_Mapped_to_Genome",
                            "Reads_Mapped_Confidently_to_Genome",
                            "Reads_Mapped_Confidently_to_Intergenic_Regions",
                            "Reads_Mapped_Confidently_to_Intronic_Regions",
                            "Reads_Mapped_Confidently_to_Exonic_Regions",
                            "Reads_Mapped_Confidently_to_Transcriptome",
                            "Reads_Mapped_Antisense_to_Gene",
                            "Fraction_Reads_in_Cells",
                            "Total_Genes_Detected",
                            "Median_UMI_Counts_per_Cell")
  
  # 2.7. Add Metadata 
  row.names(statistics) <- statistics$orig.ident 
  sce <- AddMetaData(sce, statistics)
  
  # 2.8. Save the results
  return(sce)
}  


#' diagnosisPlots
#' 
#' Create a series of diagnosis graphics for Seurat objects created with read10xData
#' @param sce, Seurat Object
#' @param sample.id, ID of the sample to use for the graphics
#' @param ID.colors.f, wrapper matching colors-sample.id. Passed to scale_fill_manual
#' @param umi.tresh.low, minimum number of umi reads
#' @param umi.tresh.high, maximum number of umi reads
#' @param feature.tresh.low, minimum number of features
#' @param feature.tresh.high, maximum number of features
#' @param mito.thresh.low, minimum number of mitochondrial reads (%)
#' @param mito.thresh.high, maximum number of mitochondrial reads (%)
#' @param postprocess, logical. Are the data already processed ? If not run basicSeurat to get the umap. If yes use the umap already there.
#' @param output.path, path/to/file.jpg

diagnosisPlots <- function(ssce = NULL, sample.id = NULL, postprocess = FALSE, output.path = NULL, ID.colors.f = ID.colors, umi.tresh.low = 1250, umi.tresh.high = 12000, feature.tresh.low = 500, feature.tresh.high = 3000, mito.thresh.low = 0, mito.thresh.high = 10){
  
  meta <- ssce@meta.data
  meta$sample.id <- sample.id
  
  # UMI Reads
  pMITO <- ggplot(data = meta) +
    geom_violin(aes(x = sample.id, y = mito.perc, fill = sample.id), alpha = 0.3) +
    geom_hline(aes(yintercept = mito.thresh.low), linetype = "dashed", color = "darkred") +
    geom_hline(aes(yintercept = mito.thresh.high), linetype = "dashed", color = "darkred") +
    scale_fill_manual(values = ID.colors.f()) +
    xlab("") +
    ylab("Mitochondrial Genes (%)") +
    theme_bw() +
    ggtitle("Mitochondrial Genes (%)") +
    scale_y_continuous(breaks = seq(0, 200, 10)) +
    theme(legend.position = "none")
  
  pUMI <- ggplot(data = meta) +
    geom_violin(aes(x = sample.id, y = nUMI, fill = sample.id), alpha = 0.3) +
    xlab("") +
    ylab("Reads") +
    scale_fill_manual(values = ID.colors.f()) +
    theme_bw() +
    ggtitle("UMI Counts") +
    geom_hline(aes(yintercept = umi.tresh.low), linetype = "dashed", color = "darkred") +
    geom_hline(aes(yintercept = umi.tresh.high), linetype = "dashed", color = "darkred") +
    scale_y_continuous(breaks = seq(0,1000000, 25000)) +
    theme(legend.position = "none")
  
  pFEATURES <- ggplot(data = meta) +
    geom_violin(aes(x = sample.id, y = nFeatures, fill = sample.id), alpha = 0.3) +
    xlab("") +
    ylab("Features") +
    scale_fill_manual(values = ID.colors.f()) +
    theme_bw() +
    ggtitle("Features Count") +
    geom_hline(aes(yintercept = feature.tresh.high), linetype = "dashed", color = "darkred") +
    geom_hline(aes(yintercept = feature.tresh.low), linetype = "dashed", color = "darkred") +
    scale_y_continuous(breaks = seq(0,20000, 1000)) +
    theme(legend.position = "none")
  
  # Cell Distribution
  pCellDistrib <- meta %>%
    mutate(outlier.type = case_when(
      isDoublet == TRUE ~ "doublet",
      percent.mito < mito.thresh.low | percent.mito > mito.thresh.high ~ "outlier.mito",
      nUMI < umi.tresh.low & nUMI > umi.tresh.high ~ "outlier.umi",
      nFeature_RNA < feature.tresh.low & nFeature_RNA > feature.tresh.high ~ "outlier.features",
      TRUE ~ "Cell"
    )) %>%
    ggplot() +
    geom_bar(aes(x = "", fill = outlier.type), color = "black") +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    scale_y_continuous(breaks = seq(0, nrow(meta)+1000, 1000)) +
    theme(legend.title =  element_blank(),
          panel.border = element_blank(),
          axis.ticks.x =element_blank(), 
          axis.text.x = element_blank(), 
          axis.line.x = element_blank(), 
          axis.title.x = element_blank(),
          panel.grid = element_blank())
  
  if(postprocess == FALSE){
    
    # Saturation
    pSaturation <- meta %>%
      summarise(saturation = unique(Sequencing_Saturation),
                read_per_cell = Number_of_Reads / Estimated_Number_of_Cells) %>%
      distinct() %>%
      ggplot() +
      geom_bar(aes(x = 1, y = 25000, fill = "A"), stat = "identity", color = "black") +
      geom_bar(aes(x = 1, y = read_per_cell, fill = "B"), stat = "identity", color = "black", alpha = 0.5) +
      geom_hline(yintercept = 25000, linetype = "dashed", color = "darkgrey") +
      xlab("") +
      ylab("Estimated Seq. Depth/Cell") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      coord_flip() +
      geom_text(aes(y = 1000, x = 1, color = "white", label = paste0("Saturation: ", saturation, "\nReads per Cell: ", floor(read_per_cell))), hjust = 0) +
      theme(legend.position = "none", 
            axis.ticks.y =element_blank(), 
            panel.border = element_blank(),
            axis.text.y = element_blank(), 
            axis.line.y = element_blank(), 
            panel.grid = element_blank())
    
    #UMAPs
    umap.ssce <- basicSeurat(ssce, SCT = TRUE, npcs = 25, regress.mito = TRUE)
    umap.filtered <- basicSeurat(ssce[,ssce$nUMI > umi.tresh.low & ssce$nUMI < umi.tresh.high & 
                                        ssce$nFeature_RNA > feature.tresh.low & ssce$nFeature_RNA < feature.tresh.high & 
                                        ssce$percent.mito > mito.thresh.low & ssce$percent.mito < mito.thresh.high], SCT = TRUE, npcs = 20, regress.mito = TRUE)
    umap.ssce$bad <- !colnames(umap.ssce) %in% colnames(umap.filtered)
    
    pUMAP <- plot_grid(nrow = 3,
                       DimPlot(umap.ssce, label = T) + theme(legend.position = "none") + ggtitle("Quick UMAP"), 
                       DimPlot(umap.ssce, group.by = "bad", cols = RColorBrewer::brewer.pal("Set1", n = 3)) + ggtitle("Quick Bad Cells"),
                       FeaturePlot(umap.ssce, features = "percent.mito") + scale_color_viridis_c(option = "viridis") +  ggtitle("Quick UMAP - Mitochondria (%)"), 
                       FeaturePlot(umap.ssce, features = "nUMI") + scale_color_viridis_c(option = "plasma") +  ggtitle("Quick UMAP - nUMI"), 
                       FeaturePlot(umap.ssce, features = "nFeature_RNA") + scale_color_viridis_c(option = "plasma") +  ggtitle("Quick UMAP - nFeatures"), 
                       DimPlot(umap.filtered, label = T) + theme(legend.position = "none") + ggtitle("Quick UMAP - Filtered"),
                       pCellDistrib)
    
    
    pdiag <- plot_grid(ncol = 1, rel_heights = c(0.2, 0.8),
                       plot_grid(pMITO, pUMI, pFEATURES, pSaturation, nrow = 1),
                       pUMAP
    ) %>% return()
    
  } else {
    
    pUMAP <- plot_grid(nrow = 3,
                       DimPlot(ssce, label = T) + theme(legend.position = "none"), 
                       DimPlot(ssce, group.by = "newID", cols = ID.colors.f()) + ggtitle("IDs"),
                       FeaturePlot(ssce, features = "percent.mito") + scale_color_viridis_c(option = "viridis") +  ggtitle("Mitochondria (%)"), 
                       FeaturePlot(ssce, features = "nUMI") + scale_color_viridis_c(option = "plasma") +  ggtitle("nUMI"), 
                       FeaturePlot(ssce, features = "nFeature_RNA") + scale_color_viridis_c(option = "plasma") +  ggtitle("nFeatures"), 
                       pCellDistrib)
    
    pdiag <- plot_grid(ncol = 1, rel_heights = c(0.2, 0.8),
                       plot_grid(pMITO, pUMI, pFEATURES, nrow = 1),
                       pUMAP
    ) 
    
  }
  
  ggsave(pdiag, filename = output.path, height = 20, width = 20, dpi = 250)
  
  return(umap.filtered)
}



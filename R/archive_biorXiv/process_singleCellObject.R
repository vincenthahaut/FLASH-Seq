library(scRNAseq)
library(scater)
library(scran)
library(Seurat) # For visualization

# 1. Load the singleCellObject
# This object contains the FS-21c and FS-LA-16c hPBMCs. 
# A prefilter was applied to remove cells <70% uniquely mapped reads and <100000 raw reads
# Two metadata columns were added (PLATE and GROUP [FS-21c or LA-16c])
mysce <- read_rds("FS_LA_sce.rds")

# 2. Remove Lowly expressed genes (=<10 cells)
mysce <- mysce[rowSums(assay(mysce)) > 10,] 

# 3. Quality control.
is.mito <- grepl("^MT-", rownames(mysce))
mysce <- addPerCellQC(mysce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(mysce, percent_subsets="subsets_Mito_percent")
mysce <- mysce[, !filtered$discard]

# 4. Remove Genes without direct interest for this study (Mitochondria, ribosomal, MALAT1)
# plotHighestExprs(mysce, exprs_values = "counts")
mysce <- mysce[!grepl("MT-|MALAT1|RPS|RPL", row.names(mysce)),] 

# 5. Library Normalization
lib.sf <- librarySizeFactors(mysce)
clust.sce <- quickCluster(mysce) 
# deconv.sf <- calculateSumFactors(mysce, cluster=clust.sce)
# hist(log10(deconv.sf), xlab="Log10[Size factor]", col='grey80')
mysce <- computeSumFactors(mysce, cluster=clust.sce, min.mean=0.1)

# 7. Log-Normalization
mysce <- logNormCounts(mysce)

# 8. Feature selection.
dec <- modelGeneVar(mysce)

# 9. Visualize the fit
fit.pbmc <- metadata(dec)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# 10. Highly variable genes (10%)
hvg <- getTopHVGs(dec, prop=0.1)

# 11. Dimensionality reduction.
# mysce <- runPCA(mysce, ncomponents=40, subset_row=hvg)
# percent.var <- attr(reducedDim(mysce), "percentVar")
# plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
set.seed(42)
mysce <- runPCA(mysce, ncomponents=10, subset_row=hvg)
mysce <- runUMAP(mysce, dimred = 'PCA', external_neighbors=TRUE)

# 12. Clustering.
g <- buildSNNGraph(mysce, k = 6, use.dimred = 'PCA')
mysce$CLUSTER <- factor(igraph::cluster_louvain(g)$membership)

# 13. Graphics
plotUMAP(mysce, colour_by="GROUP")
plotUMAP(mysce, colour_by="CLUSTER")

# 14. Cell types
mysce.seu  <- as.Seurat(mysce)
# CD45 Cells: FeaturePlot(mysce.seu, c("PTPRC"))
mysce$CELLTYPE <- case_when(
  # FeaturePlot(mysce.seu, features=c("CCR7", "SELL", "CD28"), cols = RColorBrewer::brewer.pal("Spectral", n = 8))
  mysce$CLUSTER == 3 ~ "Naïve T-Cells",
  # FeaturePlot(mysce.seu, features=c("CD4", "IL7R", "CD3D"))
  mysce$CLUSTER == 1 ~ "CD4 T-Cells",
  # FeaturePlot(mysce.seu, features=c("CD8A", "CD8B", "PRF1", "GZMA", "CD3D"))
  mysce$CLUSTER == 7 ~ "CD8+ T-Cells",
  # FeaturePlot(mysce.seu, c("FCGR3A", "NKG7", "GNLY", "NCAM1"))
  # FeaturePlot(mysce.seu, c("TRDV2", "TRGV9", "KLRD1", "KLRG1"))
  mysce$CLUSTER == 6 ~ "NK-Cells",
  # FeaturePlot(mysce.seu, c("CD79A", "JCHAIN", "MS4A1", "IGHM", "IGHD"))  
  mysce$CLUSTER == 2 ~ "B-Cells",
  # Non-Classical Monocytes
  # FeaturePlot(mysce.seu, c("CD14", "FCGR3A", "CDKN1C", "LYZ"))
  mysce$CLUSTER == 4 ~ "Non-Classical\nMonocytes",
  mysce$CLUSTER == 5 ~ "CD14+ Monocytes")

# 15. Pheatmap
rdims <- as.data.frame(reducedDim(mysce, type = "UMAP")) %>%
  rownames_to_column("ID") %>%
  left_join(select(as.data.frame(mysce@colData) %>% rownames_to_column("ID"), CLUSTER, GROUP, CELLTYPE, ID), by = "ID")
colnames(rdims) <- c("ID", "UMAP_1", "UMAP_2", "CLUSTER", "GROUP", "CELLTYPE")

annots <- data.frame(GROUP = mysce$GROUP, CELLTYPE = mysce$CELLTYPE)
row.names(annots) <- colnames(mysce)

tiff("PBMC_PHEATMAP.tiff", width = 1800, height = 900)
plotHeatmap(
  mysce,
  features = c("CD14",  "LYZ", "CDKN1C", "CD79A", "FCER1A", "CD1C",
               "JCHAIN", "MS4A1", "IGHM", "IGHD",
               "NKG7", "GNLY", "NCAM1", "FCGR3A", "PRF1", "GZMA", 
               "TRDV2", "TRGV9", "TRDC", 
               "CD8A", "CD8B", "CD3E",
               "TCF7", "CD4", "IL7R", "CCR7", "LRRN3", "LEF1", "NOG", "PRDM1"), 
  annotation = annots,
  fontsize = 18,  
  colour_columns_by = c("GROUP", "CELLTYPE")
)
dev.off()

# 16. UMAP (groups and clusters)
p.groups <- ggplot(rdims, aes(x = UMAP_1, y = UMAP_2, color = GROUP, fill = GROUP)) +
  geom_point(alpha = 0.5, size = 1.25) +
  theme_cowplot(font_size = 20) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#B2435F", "#43b296"))

p.celltypes <- ggplot(rdims, aes(x = UMAP_1, y = UMAP_2, color = CELLTYPE)) +
  geom_point(alpha = 0.5, size = 1.2) +
  theme_cowplot(font_size = 20) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = mycolors) +
  guides(color = guide_legend(override.aes = list(size=3)))

ggsave(p.celltypes + theme(legend.position = "none"), filename = "PBMC_cellTypes.tiff", width = 6, height = 4, dpi = 300)
ggsave(p.groups + theme(legend.text = element_text(size = 14)), filename = "PBMC_groups.tiff", width = 6, height = 4, dpi = 300)




#!/bin/R

# 0. Functions 
basicSeurat <- function(mysce = NULL, select.ids = mydata$sce.ids, GROUP = NULL, ndims = 15, nfeatures = 3000, out.path = "/home/vincent.hahaut/PBMC_ssce.rds"){
    library(Seurat)
    ssce <- Seurat::CreateSeuratObject(counts = as.matrix(assay(mysce[,colnames(mysce) %in% select.ids])))
    # ssce <- ssce[rowSums(Seurat::GetAssayData(ssce)) > 5,] # Remove lowly expressed genes
    ssce$PLATE <- str_replace_all(str_split(colnames(ssce), "_", simplify = T)[,1], pattern = "X|_", "")
    ssce$GROUP <- GROUP
    ssce[["percent.mt"]] <- PercentageFeatureSet(ssce, pattern = "^MT-")
    ssce <- ssce[!grepl(row.names(ssce), pattern = "MT-|MALAT1|RPS|RPL"),]
    ssce <- Seurat::SCTransform(ssce, vars.to.regress = c("nCount_RNA", "percent.mt"), variable.features.n = nfeatures, seed.use = 42)
    ssce <- Seurat::RunPCA(ssce, npcs = 40)
    ssce <- Seurat::RunUMAP(ssce, dims = 1:ndims, umap.method = "uwot")
    ssce <- Seurat::FindNeighbors(ssce, dims = 1:ndims)
    ssce <- Seurat::FindClusters(ssce, resolution = 1.5)
    write_rds(ssce, out.path)
    return(ssce)
}


# 1. Load the count matrix and select the cells to process 
sce <- read_rds("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.2021-12-05.rds.bz")
ID.SS3 <- grepl(colnames(sce), pattern = "SS3")
ID.FS <- grepl(colnames(sce), pattern = "FS")

# 2. Process cells
ssce.ss3 <- basicSeurat(mysce = sce, GROUP = "SS3", select.ids = ID.SS3, ndims = 8, nfeatures = 3000, out.path = "/home/vincent.hahaut/PBMC_SS3_ssce.rds")
ssce.only <- basicSeurat(mysce = sce, GROUP = "FS - 25 µL", select.ids = ID.FS, ndims = 12, nfeatures = 3000, out.path = "/home/vincent.hahaut/PBMC_FS_ssce.rds")

# 3. Add Azimuth predictions
# Generated online using the PBMC_SS3_ssce.rds / PBMC_FS_ssce.rds objects
ss3.ex.pred <- read_tsv("/home/vincent.hahaut/Downloads/SS3_ex_azimuth.tsv")
FS.ex.pred <- read_tsv("/home/vincent.hahaut/Downloads/FS_ex_azimuth.tsv")

predictions <- bind_rows(ss3.ex.pred,
          FS.ex.pred) %>% as.data.frame()

row.names(predictions) <- predictions$cell
predictions$cell <- NULL
ssce.ss3 <- Seurat::AddMetaData(ssce.ss3, predictions)
ssce.only <- Seurat::AddMetaData(ssce.only, predictions)
ssce.ss3$predicted.celltype.l2 <- str_replace_all(ssce.ss3$predicted.celltype.l2, " ", "_")
ssce.only$predicted.celltype.l2 <- str_replace_all(ssce.only$predicted.celltype.l2, " ", "_")

# 4. Explore the results
predicted.l2.colors <- c(
  "B_intermediate" = "#6B8993FF",
  "B_memory" = "#D95F02",
  "B_naive" = "#7570B3",
  "Plasmablast" = "#E7298A",
  "CD4_CTL" = "#9C7E70FF",
  "CD4_Naive" = "#E6AB02",
  "CD4_Proliferating" = "#A6761D",
  "CD4_TCM" = "#AEADB0FF",
  "CD4_TEM" = "#0C5BB0FF",
  "Treg" = "#94D4D4FF",
  "CD8_Naive" = "#15983DFF",
  "CD8_Proliferating" = "#EC579AFF",
  "CD8_TCM" = "#A1C720FF",
  "CD8_TEM" = "black",
  "cDC2" = "#FA6B09FF", 
  "pDC" = "#9DDAF5FF",
  "CD14_Mono" = "#66A61E", 
  "CD16_Mono" = "#9A703EFF",
  "NK" = "#EE0011FF", 
  "Eryth" = "#666666",
  "Platelet" = "#16A08CFF",
  "gdT" = "#9966CC",
  "MAIT" = "#11776CFF")

ggsave(plot_grid(
  DimPlot(ssce.ss3, group.by = "predicted.celltype.l2", cols = mycolors) + 
    theme(plot.title = element_blank(), legend.position = "none") + 
    scale_color_manual(values = predicted.l2.colors),
  DimPlot(ssce.only, group.by = "predicted.celltype.l2", cols = mycolors) + 
    theme(plot.title = element_blank(), legend.position = "none") + 
    scale_color_manual(values = predicted.l2.colors)),
bg = "white", filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/DimPlot_PBMC_FS.predicted.tiff", dpi = 300, width = 12, height = 6)

# 5. Check up:
B_intermediate = c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B", "CD3D")
B_memory =c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781", "")
B_naive =c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3", "")
Plasmablast =c("IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","TNFRSF13B","POU2AF1","CPNE5","HRASLS2","NT5DC2", "")
CD4_CTL =c("GZMH","CD4","FGFBP2","ITGB1","GZMA","CST7","GNLY","B2M","IL32","NKG7", "")
CD4_Naive =c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1", "CD3D")
CD4_Proliferating =c("MKI67","TOP2A","PCLAF","CENPF","TYMS","NUSAP1","ASPM","PTTG1","TPX2","RRM2", "")
CD4_TCM =c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL", "")
CD4_TEM =c("IL7R","CCL5","FYB1","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3", "")
Treg =c("RTKN2","FOXP3","AC133644.2","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2", "")
CD8_Naive =c("CD8B","S100B","CCR7","RGS10","NOSIP","LINC02446","LEF1","CRTAM","CD8A","OXNAD1", "CD3D")
CD8_Proliferating =c("MKI67","CD8B","TYMS","TRAC","PCLAF","CD3D","CLSPN","CD3G","TK1","RRM2", "")
CD8_TCM =c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB", "")
CD8_TEM =c("CCL5","GZMH","CD8A","TRAC","KLRD1","NKG7","GZMK","CST7","CD8B","TRGC2", "")
# ASDC =c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","DNASE1L3","CLEC4C","GAS6")
# cDC1 =c("CLEC9A","DNASE1L3","C1orf54","IDO1","CLNK","CADM1","FLT3","ENPP1","XCR1","NDRG2")
cDC2 =c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN","SLC38A1","NDRG2","AFF3", "")
pDC =c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3", "")
CD14_Mono =c("S100A9","CTSS","S100A8","LYZ","VCAN","S100A12","IL1B","CD14","G0S2","FCN1", "")
CD16_Mono =c("CDKN1C","FCGR3A","PTPRC","LST1","IER5","MS4A7","RHOC","IFITM3","AIF1","HES4", "")
NK =c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1", "CD3D")
# NK_Proliferating =c("MKI67","KLRF1","TYMS","TRDC","TOP2A","FCER1G","PCLAF","CD247","CLSPN","ASPM")
# NK_CD56bright =c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A")
Eryth =c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","IFIT1B","TRIM58","SELENBP1","TMCC2", "")
# HSPC =c("SPINK2","PRSS57","CYTL1","EGFL7","GATA2","CD34","SMIM24","AVP","MYB","LAPTM4B"),
Platelet =c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9", "")
# ILC =c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS"),
#dnT =c("PTPN3","MIR4422HG","NUCB2","CAV1","DTHD1","GZMA","MYB","FXYD2","GZMK","AC004585.1"),
gdT =c("TRDC","TRGC1","TRGC2","KLRC1","NKG7","TRDV2","CD7","TRGV9","KLRD1","KLRG1", "CD3D")
MAIT =c("KLRB1","NKG7","GZMK","IL7R","SLC4A10","GZMA","CXCR6","PRSS35","RBM24","NCR3", "")

mysce <- ssce.only
FeaturePlot(mysce, B_naive)
FeaturePlot(mysce, Plasmablast)
FeaturePlot(mysce, B_intermediate)
FeaturePlot(mysce, CD4_Naive)
FeaturePlot(mysce, CD4_Proliferating)
FeaturePlot(mysce, CD4_TCM)
FeaturePlot(mysce, CD4_TEM)

FeaturePlot(mysce, CD8_Naive)
FeaturePlot(mysce, CD8_Proliferating)
FeaturePlot(mysce, CD8_TCM)
FeaturePlot(mysce, CD8_TEM)
FeaturePlot(mysce, gdT)
FeaturePlot(mysce, MAIT)
FeaturePlot(mysce, NK)

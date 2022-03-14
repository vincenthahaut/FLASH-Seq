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
ID.FS <- grepl(colnames(sce), pattern = "FS")

# 2. Process cells
ssce.la <- basicSeurat(mysce = sce, GROUP = "FS-LA", select.ids = ID.FSLA, ndims = 12, nfeatures = 3000, out.path = "/home/vincent.hahaut/Desktop/BMC_LA_ssce.rds")

# 2. Add Groups
ssce.la$GROUP <- ifelse(grepl("X255|X241|X287|X290|X291", colnames(ssce.la)), "FS-LA", "FS")
ssce.la$PLATE <- ssce.la$orig.ident

# 3. Explore Results
Seurat::DimPlot(ssce.la, group.by = "GROUP")
Seurat::DimPlot(ssce.la, group.by = "PLATE")
Seurat::DimPlot(ssce.la, group.by = "seurat_clusters")

# 4. Azimuth Annotation
# Generated online using the object "/home/vincent.hahaut/Desktop/PBMC_LA_ssce.rds"
pred <- as.data.frame(read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.tsv"))
row.names(pred) <- pred$cell
ssce.la <- AddMetaData(ssce.la, pred)

# 5. Remove some labels 
# Too few cells and too unsure
ssce.la$predicted.celltype.l2 <- ifelse(ssce.la$predicted.celltype.l2 %in% c("dnT", "HSPC", "NK CD56bright", "NK Proliferating", "CD4 CTL"), NA, ssce.la$predicted.celltype.l2)
ssce.la$predicted.celltype.l2 <- str_replace_all(ssce.la$predicted.celltype.l2, " ", "_")

# 6. Explore Results 2
predicted.l2.colors <- c("B_intermediate" = "#6B8993FF",
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

plot_grid(Seurat::DimPlot(ssce.la, group.by = "predicted.celltype.l2", label = T, cols = mycolors) + scale_color_manual(values = predicted.l2.colors),
          Seurat::DimPlot(ssce.la, group.by = "GROUP"),
          Seurat::DimPlot(ssce.la, group.by = "PLATE"),
          Seurat::DimPlot(ssce.la, group.by = "seurat_clusters", label = T))

# markers <- FindAllMarkers(ssce.la, logfc.threshold = 0.75, min.pct = 0.2, test.use = "MAST")

# 7. Save Graphics
ggsave(filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/DimPlot_FS_LA.tiff", 
       plot_grid(rel_widths = c(0.4,0.4, 0.2),
  Seurat::DimPlot(ssce.la, group.by = "GROUP", pt.size = 1.25, cols = c("#FF850D", "#7f3fb2")),
  Seurat::DimPlot(ssce.la, group.by = "predicted.celltype.l2", pt.size = 1.25, cols = mycolors) + 
    scale_color_manual(values = predicted.l2.colors) +
    theme(legend.position = "none"),
  get_legend(Seurat::DimPlot(ssce.la, group.by = "predicted.celltype.l2", pt.size = 1.25, cols = mycolors) + 
    scale_color_manual(values = predicted.l2.colors)), nrow = 1), dpi = 300, width = 16, height = 6, bg = "white")

ggsave(filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/DimPlot_FS_LA.predicted.tiff", bg = "white",  width = 8, height = 8,
  Seurat::DimPlot(ssce.la, group.by = "predicted.celltype.l2", pt.size = 1.25, cols = mycolors) + 
    scale_color_manual(values = predicted.l2.colors) +
    theme(legend.position = "none"))

# 8. Heatmap
ssce.la <- ScaleData(ssce.la, features = row.names(ssce.la))
makeHeatmap <- function(mysce = NULL, size = 5, cell_type = NULL){
  markers <- data_frame(
    B_intermediate = c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B", "CD3D"),
    B_memory =c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781", ""),
    B_naive =c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3", ""),
    Plasmablast =c("IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","TNFRSF13B","POU2AF1","CPNE5","HRASLS2","NT5DC2", ""),
    CD4_CTL =c("GZMH","CD4","FGFBP2","ITGB1","GZMA","CST7","GNLY","B2M","IL32","NKG7", ""),
    CD4_Naive =c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1", "CD3D"),
    CD4_Proliferating =c("MKI67","TOP2A","PCLAF","CENPF","TYMS","NUSAP1","ASPM","PTTG1","TPX2","RRM2", ""),
    CD4_TCM =c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL", ""),
    CD4_TEM =c("IL7R","CCL5","FYB1","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3", ""),
    Treg =c("RTKN2","FOXP3","AC133644.2","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2", ""),
    CD8_Naive =c("CD8B","S100B","CCR7","RGS10","NOSIP","LINC02446","LEF1","CRTAM","CD8A","OXNAD1", "CD3D"),
    CD8_Proliferating =c("MKI67","CD8B","TYMS","TRAC","PCLAF","CD3D","CLSPN","CD3G","TK1","RRM2", ""),
    CD8_TCM =c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB", ""),
    CD8_TEM =c("CCL5","GZMH","CD8A","TRAC","KLRD1","NKG7","GZMK","CST7","CD8B","TRGC2", ""),
    # ASDC =c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","DNASE1L3","CLEC4C","GAS6"),
    # cDC1 =c("CLEC9A","DNASE1L3","C1orf54","IDO1","CLNK","CADM1","FLT3","ENPP1","XCR1","NDRG2"),
    cDC2 =c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN","SLC38A1","NDRG2","AFF3", ""),
    pDC =c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3", ""),
    CD14_Mono =c("S100A9","CTSS","S100A8","LYZ","VCAN","S100A12","IL1B","CD14","G0S2","FCN1", ""),
    CD16_Mono =c("CDKN1C","FCGR3A","PTPRC","LST1","IER5","MS4A7","RHOC","IFITM3","AIF1","HES4", ""),
    NK =c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1", "CD3D"),
    # NK_Proliferating =c("MKI67","KLRF1","TYMS","TRDC","TOP2A","FCER1G","PCLAF","CD247","CLSPN","ASPM"),
    # NK_CD56bright =c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A"),
    Eryth =c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","IFIT1B","TRIM58","SELENBP1","TMCC2", ""),
    # HSPC =c("SPINK2","PRSS57","CYTL1","EGFL7","GATA2","CD34","SMIM24","AVP","MYB","LAPTM4B"),
    Platelet =c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9", ""),
    # ILC =c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS"),
    #dnT =c("PTPN3","MIR4422HG","NUCB2","CAV1","DTHD1","GZMA","MYB","FXYD2","GZMK","AC004585.1"),
    gdT =c("TRDC","TRGC1","TRGC2","KLRC1","NKG7","TRDV2","CD7","TRGV9","KLRD1","KLRG1", "CD3D"),
    MAIT =c("KLRB1","NKG7","GZMK","IL7R","SLC4A10","GZMA","CXCR6","PRSS35","RBM24","NCR3", "")
  )

  DoHeatmap(mysce, 
          cells = colnames(mysce)[mysce$predicted.celltype.l2 %in% cell_type],
          features = as.vector(unlist(markers[cell_type])), 
          group.by = "predicted.celltype.l2", group.bar = T, raster = T, angle = 90, size = size, group.colors = predicted.l2.colors) +
    guides(color=FALSE) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(name = "RdBu", 8))) %>%
    return()

} 

# 9. Individual heatmaps
# Not included in the final manuscript
p.heat.B <- makeHeatmap(mysce = ssce.la, cell_type = c("B_intermediate", "B_memory", "B_naive"), size = 0)
p.heat.CD4 <- makeHeatmap(mysce = ssce.la, cell_type = c("CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg"), size = 0)
p.heat.CD8 <- makeHeatmap(mysce = ssce.la, cell_type = c("CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM"), size = 0)
p.heat.DC <- makeHeatmap(mysce = ssce.la, cell_type = c("cDC2", "pDC"), size = 0)
p.heat.mono <- makeHeatmap(mysce = ssce.la, cell_type = c("CD14_Mono", "CD16_Mono"), size = 0)
p.heat.nk <- makeHeatmap(mysce = ssce.la, cell_type = c("NK"), size = 0)
p.heat.otherT <- makeHeatmap(mysce = ssce.la, cell_type = c("MAIT", "gdT"), size = 0)

p.heat.all <- makeHeatmap(mysce = ssce.la, cell_type = 
                               c("B_intermediate", "B_memory", "B_naive","CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg","MAIT", "gdT","CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM","CD14_Mono", "CD16_Mono","NK","cDC2", "pDC"))

ggsave(p.heat.all, bg = "white", width = 18, height = 16, filename =  "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/FS_LA_heatmap_all.tiff", dpi = 300)
ggsave(plot_grid(p.heat.B, p.heat.CD4, p.heat.CD8, p.heat.DC, 
                 p.heat.mono, p.heat.nk, p.heat.otherT, ncol = 2), bg = "white", width = 12, height = 18, filename =  "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/FS_LA_heatmap_sep.tiff", dpi = 300)


# 10. Label statibility
# For each downsampling, run a separated analysis and then use azimuth to get the cell types
ssce.la.5K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_5000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.5K.rds")
ssce.la.10K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_10000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.10K.rds")
ssce.la.20K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_20000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.20K.rds")
ssce.la.40K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_40000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.40K.rds")
ssce.la.50K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_50000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.50K.rds")
ssce.la.75K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_75000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.75K.rds")
ssce.la.100K <- basicSeurat(mysce = sce, GROUP = "ALL", select.ids = str_replace(mydata$sce.ids, "_125000_1", "_100000_1"), ndims = 10, nfeatures = 2000, out.path = "/home/vincent.hahaut/PBMC_LA_ssce.100K.rds")

pred.5k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.5k.tsv") %>% 
  mutate(DEPTH = "5K", 
         ID = str_replace(cell, "_5000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.10k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.10k.tsv") %>% 
  mutate(DEPTH = "10K", 
         ID = str_replace(cell, "_10000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.20k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.20k.tsv") %>% 
  mutate(DEPTH = "20K", 
         ID = str_replace(cell, "_20000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.40k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.40k.tsv") %>% 
  mutate(DEPTH = "40K", 
         ID = str_replace(cell, "_40000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.50k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.50k.tsv") %>% 
  mutate(DEPTH = "50K", 
         ID = str_replace(cell, "_50000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.75k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.75k.tsv") %>% 
  mutate(DEPTH = "75K", 
         ID = str_replace(cell, "_75000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.100k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.100k.tsv") %>% 
  mutate(DEPTH = "100K", 
         ID = str_replace(cell, "_100000_1", "")) %>%
  select(ID, predicted.celltype.l2, DEPTH)
pred.125k <- read_tsv("/home/vincent.hahaut/Downloads/FS_LA_azimuth.tsv") %>% 
  mutate(DEPTH = "125K", 
         ID = str_replace(cell, "_125000_1", "")) %>%
  select(ID, predicted.celltype.l2)

p.label.stability <- bind_rows(pred.5k, pred.10k, pred.20k, pred.40k, pred.50k, pred.75k, pred.100k) %>%
  left_join(pred.125k, by = "ID") %>%
  mutate(overlap = predicted.celltype.l2.x == predicted.celltype.l2.y,
         DEPTH = factor(DEPTH, levels = c("5K", "10K", "20K", "40K", "50K", "75K", "100K", "125K")),
         overlap = ifelse(overlap == TRUE, "Matches Reference Cell Type", "Diverges from Reference Cell Type"),
         GROUP = ifelse(grepl(ID, pattern = "X255|X241|X287|X290|X291"), "FS-LA", "FS" )) %>%
  group_by(DEPTH, GROUP) %>%
  mutate(total = n()) %>%
  group_by(DEPTH, overlap, GROUP) %>%
  summarise(perc = 100*n()/unique(total)) %>%
  ggplot(aes(x = DEPTH, y = perc, fill = overlap)) + 
  geom_bar(stat="identity", color = "black") +
  geom_hline(yintercept = 91, linetype = "dashed", color = "lightgrey") +
  scale_y_continuous(breaks = c(91, seq(0,85,10), 100)) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_fill_brewer(palette = "Set1") +
  xlab("Downsampled Reads") + 
  ylab("Cell (%)") + 
  facet_wrap("GROUP")

ggsave(p.label.stability, bg = "white", width = 8, height = 6, filename =  "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/FS_LA_label.statiblity.tiff", dpi = 300)


# 11. TRUST4

# 11.0. function to parse TRUST4 results
# Get the best assembly
# Select only "TRBC", "TRAC",  "TRGC",  "TRDC" chains
# Remove the allele info
# tidy the results
parseTrust4 <- function(x = NULL, DEPTH = NULL){
  x %>%
    mutate(cid = str_replace(pattern = "_assemble.+$", replacement = "", cid)) %>%
    filter(CDR3aa != "out_of_frame") %>%
    group_by(cid) %>%
    arrange(desc(`#count`)) %>%
    filter(row_number()==1) %>%
    mutate(reads = `#count`) %>%
    dplyr::select(CDR3aa, V, D, J, C, cid, cid_full_length, reads) %>%
    filter(C %in% c("TRBC", "TRAC",  "TRGC",  "TRDC")) %>%
    ungroup() %>% 
    mutate(V = str_replace(V, "\\*.+$", ""),
           D = str_replace(D, "\\*.+$", ""),
           J = str_replace(J, "\\*.+$", ""),
           DEPTH = DEPTH,
           TCR = ifelse(cid_full_length == 1, paste0(V, "_", J, "_", C, "_", CDR3aa), "")) %>%
    return()
}

# 11.1 Read Trust4 results
trust.5K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_5K_report.tsv") %>% parseTrust4(x = ., DEPTH = 5000)
trust.10K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_10K_report.tsv") %>% parseTrust4(x = ., DEPTH = 10000)
trust.20K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_20K_report.tsv") %>% parseTrust4(x = ., DEPTH = 20000)
trust.40K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_40K_report.tsv") %>% parseTrust4(x = ., DEPTH = 40000)
trust.50K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_50K_report.tsv") %>% parseTrust4(x = ., DEPTH = 50000)
trust.75K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_75K_report.tsv") %>% parseTrust4(x = ., DEPTH = 75000)
trust.100K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_100K_report.tsv") %>% parseTrust4(x = ., DEPTH = 100000)
trust.125K <- read_tsv("/home/vincent.hahaut/data_storage/TRUST4/TRUST_125K_report.tsv") %>% parseTrust4(x = ., DEPTH = 125000)

trust4 <- list(trust.5K, trust.10K, trust.20K, trust.40K, trust.50K, trust.75K, trust.100K, trust.125K,)

# 11.2. Aggregate results
# Compare the TCR detection between 125K and the rest
ref <- dplyr::select(trust.125K, cid, TCR, CDR3aa) %>%
  rename("TCR" = "REF.TCR", "CDR3aa" = "REF.CDR3aa")
  
trust4.aggr <- bind_rows(trust4) %>%
  left_join(ref) %>%
  filter(nchar(CDR3aa) >= 8, !is.na(CDR3aa), CDR3aa != "out_of_frame") %>%
  filter(!is.na(REF.TCR)) %>%
  mutate(`V(D)JC-CDR3aa TCR` = TCR == REF.TCR, 
         `CDR3aa Only` = CDR3aa == REF.CDR3aa) %>%
  dplyr::select(DEPTH, TCR.match, CDR3aa.match) %>%
  pivot_longer(-DEPTH) %>%
  mutate(value = ifelse(value == TRUE, value, FALSE),
         FILL = "TCR")

p.trust4.full.tcr <- ggplot(trust4.aggr, aes(x = as.factor(DEPTH), fill = value)) + 
  geom_hline(yintercept = seq(0,650,50), linetype = "dashed", color = "darkgrey") +
  geom_bar(color = "black") +
  facet_wrap("name", nrow = 2) +
  ylab("Cells") +
  xlab("Downsampled Reads") + 
  theme_bw(base_size = 20) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(labels = names(trust4)) +
  scale_fill_brewer(palette = "Set1") + 
  ggtitle("Valid TCR - overlaping with REF")

ggsave(p.trust4.full.tcr, filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/TRUST4_TCRcounts.tiff", dpi = 300, height = 8, width = 9, bg = "white")


# 11.3. Visualize TCR distribution on UMAP
ref.tcr <- trust.125K %>% 
  group_by(TCR) %>%
  mutate(n = sum(nchar(TCR) > 10)) %>% 
  as.data.frame() %>%
  mutate(ID.tcr = str_replace_all(paste0("X", trust.125K$cid, "_125000_1"), "-", "\\."),
         ID.tcr = ifelse(grepl("X292", ID.tcr), str_replace(ID.tcr, "_S.+_R1_001", ""), ID.tcr),
         is.TCR = TRUE)
row.names(ref.tcr) <- ref.tcr$ID.tcr

ssce.la <- AddMetaData(ssce.la, ref.tcr)
ssce.la$is.full.TCR<- nchar(ssce.la$TCR) > 10
ssce.la$is.expanded <- ssce.la$n > 3
ssce.la$is.CDR3 <- nchar(ssce.la$CDR3aa) > 8 & ssce.la$CDR3aa != "out_of_frame"

dimplot <- plot_grid(
  DimPlot(ssce.la, group.by = "is.TCR", pt.size = 1.25, cols = c("darkred")),
  DimPlot(ssce.la, group.by = "is.full.TCR", pt.size = 1.25, cols = c("grey", "darkblue")), nrow = 1)

ggsave(plot_grid(p.trust4, dimplot, nrow = 2), filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/TRUST4.tiff", dpi = 300, height = 12, width = 16, bg = "white")

sum(ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T") & ssce.la$is.TCR, na.rm = T) / sum(ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T"), na.rm = T)
sum(!ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T") & ssce.la$is.TCR, na.rm = T)
sum(ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T"), na.rm = T)
sum(ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T") & ssce.la$is.CDR3, na.rm = T) / sum(ssce.la$predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "other T"), na.rm = T)

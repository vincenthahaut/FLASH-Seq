# 1. Load the SingleCellObject
# This object contains the gene expression for each cell
sce <- read_rds("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.2021-06-19.rds.bz")

# 2. Select the cells of interest
mysce <- sce[,colnames(sce) %in% cellIDs]

# 3. Add the GROUP info (SS2, SS3 or FS)
mysce$GROUP <- GROUP

# 4. For each group, get the expressed genes
expressedSet.genes <- sapply(c("SS2", "SS3", "FS"), function(x) 
  row.names(mysce)[rowSums(counts(mysce)[,which(mysce$GROUP == x)]) > 0])

expressedSet.genes <- intersect(intersect(expressedSet.genes$SS2, expressedSet.genes$SS3), expressedSet.genes$FS)

# 5. Subset the singleCellObject
sce.expressed <- mysce[row.names(mysce) %in% expressedSet.genes,]

# 6. Cell-to-Cell Correlations
corrs <- data_frame()
for(j in unique(c("SS2", "SS3", "FS"))){
  sce.tmp <- counts(sce.expressed[,sce.expressed$GROUP == j])
  p <- pcaPP::cor.fk(sce.tmp) 
  
  corrs <- bind_rows(corrs,
                     data_frame(GROUP = j,
                                tau = p[upper.tri(p)]))
}

# 7. Graphics
p6 <- corrs %>%
  ggplot(aes(x = reorder(GROUP, tau, median), y = tau)) + 
  geom_hline(yintercept = seq(0.55, 0.8, 0.05), linetype = "dotted", color = "darkgrey") +
  geom_boxplot(aes(fill = GROUP), width = 0.6) +
  theme_cowplot(font_size = 28) +
  theme(legend.position = "none",
        axis.text = element_text(size = 22), 
        axis.title = element_text(size = 24)) +
  xlab("") +
  ylab("Cell-to-Cell Correlation (tau)") 

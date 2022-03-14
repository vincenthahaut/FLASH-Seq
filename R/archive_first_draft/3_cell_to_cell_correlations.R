# Example of cell-to-cell correlations presented in the paper

### Cell-to-cell correlation: HEK FS, SS3, SS2

# 1. Load the SingleCellObject
# This object contains the gene expression for each cell
sce <- read_rds("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.rds.bz")

# 2. Select the cells of interest
mysce <- sce[,colnames(sce) %in% cellIDs]

# 3. Add the GROUP info (SS2, SS3 or FS)
mysce$GROUP <- GROUP

# 4. For each group, get the expressed genes
# Used a gene set common to all methods to allow comparison between the tau
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


### Cell-to-cell correlation: HEK FS-LA (various PCR cycles)

# 1. Load the SingleCellObject
# This object contains the gene expression for each cell
sce <- read_rds("/home/vincent.hahaut/data_storage/FS_intermediate_file/downsampling_featureCounts.exon.all.rds.bz")

# 2. Select the cells of interest
mysce <- sce[,colnames(sce) %in% cellIDs]

# 3. Add the GROUP info (SS2, SS3 or FS)
# "FS-19c", "FS-12c", "FS-10c", "FS-4c", "FS-8c", "FS-6c"
mysce$GROUP <- GROUP

# 4. Subset the genes expressed in > 2 cells 
expressed.genes <- lapply(unique(sce.expressed$GROUP), function(x) row.names(
  sce.expressed[,sce.expressed$GROUP == x])[rowSums(counts(sce.expressed[,sce.expressed$GROUP == x]) > 0) > 2])
expressed.genes <- Reduce(union, expressed.genes)

sce.expressed <- sce.expressed[row.names(sce.expressed) %in% expressed.genes,]

# 5. Compute pairwise group-by-group correlations
corrs <- data_frame()
for(j in unique(c("FS-19c", "FS-12c", "FS-10c", "FS-4c", "FS-8c", "FS-6c"))){
    sce.tmp <- counts(sce.expressed[,sce.expressed$GROUP == j])
    p <- pcaPP::cor.fk(sce.tmp) 
    
    corrs <- bind_rows(corrs,
                       data_frame(GROUP = j,
                                  tau = p[upper.tri(p)]))
}


# 6. Graphics
shapiro.test(corrs$tau)
corrs$GROUP <- factor(corrs$GROUP)
corrs$GROUP <- reorder(corrs$GROUP, corrs$tau, median)
du <- FSA::dunnTest(tau ~ GROUP, data=corrs, method="bonferroni")
comparisons <- str_split(du$res$Comparison, " - ", simplify = TRUE)
geomseg <- data_frame(
  x = comparisons[,1],
  xend = comparisons[,2],
  padj = padj <- du$res$P.adj,
  lab = case_when(padj < 0.05 & padj >= 0.005 ~ "*",
                  padj < 0.005 & padj >= 0.0005 ~ "**",
                  padj < 0.0005 ~ "***",
                  padj >= 0.05 ~ "ns",
                  )) %>%
  filter(lab != "ns") %>%
  mutate(y = seq(0.8, 0.8+0.01*13, 0.01)) %>%
  rowwise() %>%
  mutate(x.mid = (which(levels(corrs$GROUP) == xend) + which(levels(corrs$GROUP) == x))/2)

p6 <- corrs %>%
  ggplot(aes(x = GROUP, y = tau)) + 
  geom_hline(yintercept = seq(0.40, 0.8, 0.05), linetype = "dotted", color = "darkgrey") +
  geom_violin(aes(fill = GROUP)) +
  geom_boxplot(aes(fill = GROUP), width = 0.1) +
  theme_cowplot(font_size = 24) +
    theme(legend.position = "none",
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
        axis.title = element_text(size = 20)) +
  xlab("") +
  scale_fill_manual(values = c("FS-4c" = "#FED976", 
                               "FS-8c" = "#FEB24C", 
                               "FS-6c" = "#FD8D3C", 
                               "FS-10c" = "#FC4E2A", 
                               "FS-12c" = "#E31A1C", 
                               "FS-19c" = "#B10026")) +
  ylab("Cell-to-Cell Correlation (tau)") +
  geom_segment(data = geomseg, aes(x = x, xend = xend, y = y, yend = y)) +
  geom_text(data = geomseg, aes(x = x.mid, y = y+0.0005, label=lab), size=6)

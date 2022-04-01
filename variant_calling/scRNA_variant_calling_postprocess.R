#!/bin/R

# 0. Prerequisits
options(scipen=999)
library(tidyverse)
library(cowplot)
mycolors <- as.vector(c(yarrr::piratepal(palette = "basel"), 
                        yarrr::piratepal(palette = "info"), 
                        yarrr::piratepal(palette = "pony"), 
                        yarrr::piratepal(palette = "eternal"), 
                        yarrr::piratepal(palette = "espresso")))

gtf <- as.data.frame(rtracklayer::import("/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"))


# 1. Create a sample sheet containing:
# VCF paths
# Sample ID (ID_DEPTH)
# Unique ID (ID)
# Downsampling depth (DEPTH)
vcf.paths <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/", pattern = ".vcf", recursive = TRUE, full.names = TRUE)
vcf.sp <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/", pattern = ".vcf", recursive = TRUE, full.names = FALSE) 

pattern <- "GATK"
vcf.sp <- vcf.sp[grepl(vcf.paths, pattern = paste0("\\.", pattern, "\\."))]
vcf.paths <- vcf.paths[grepl(vcf.paths, pattern = paste0("\\.", pattern, "\\."))]

vcf.sp <- data_frame(PATHS = vcf.paths,
                     ID = str_split(vcf.sp, "/", simplify = TRUE)[,1],
                     ID_original = vcf.sp,
                     DEPTH = str_split(vcf.sp, "\\.", simplify = T)[,3])

# 2. Create a repertoire of scRNA variants shared in >10 cells
# Can be used to filter out false positives (typically -30%)
# Not used in the manuscript
vcf.sp.top <- vcf.sp %>% group_by(ID) %>% filter(DEPTH == max(as.numeric(DEPTH)))

var.dict <- list()
for(x in vcf.sp.top$PATHS){
  var.dict[[x]] <- vcfR::vcfR2tidy(vcfR::read.vcfR(x), format_fields = c("GT", "DP"))$fix %>%
    filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
    filter(DP > 2) %>%
    summarise(VARIANT = paste0(ChromKey, ":", POS, "-", REF, "/", ALT))
}

vars.dict <- bind_rows(var.dict) %>% 
  as_tibble() %>%
  group_by(VARIANT) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

ggplot(filter(vars.dict, n > 10), aes(x = as.numeric(n))) + 
  geom_histogram(bins = 150) +
  xlab("Number of Cells") + 
  theme_cowplot(font_size = 20) +
  ylab("Detected Variants") +
  scale_x_continuous(breaks = seq(10,350,25), limits = c(10, max(vars.dict$n))) 

# 3. Load the reference exome-sequencing
ref <- vcfR::read.vcfR("/home/vincent.hahaut/data_storage/ORGANOIDS/ex.hg38.vcf") %>%
  vcfR::vcfR2tidy(format_fields = c("GT", "DP"))
ref <- filter(ref$fix, nchar(REF) == 1 & nchar(ALT) == 1) %>%
  dplyr::select(ChromKey,CHROM,POS,REF,ALT,QUAL, FILTER,DP) %>%
  mutate(VARIANT = paste0(ChromKey, ":", POS, "-", REF, "/", ALT)) %>%
  filter(FILTER == "PASS")

# 4. Load the bed file containing the regions analysed in exome-seq
bed.exome <- read_tsv("/home/vincent.hahaut/Desktop/Twist_ComprehensiveExome_targets_hg38.bed", col_names = F) %>%
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "X1", start.field = "X2", end.field = "X3")
  
# 5. processVCF function
#' @param path, character, /path/to/file.vcf
#' @param ID, character, sample_id
#' @param DEPTH, numeric, sample_downsampling_depth
#' @param ref.pass, data_frame, reference exome
#' @param bed.exome, GRange object, positions explored in the ref.pass file
#' @param common.vars, data_frame, data frame containing the variants in >10 cells
processVCF <- function(paths = NULL, ID = NULL, DEPTH = NULL, ref.pass = ref, bed.exome = bed.exome, common.vars = vars.dict){
  
  # Check if the BAM file is also present
  # Will be used to determine the coverage of the exome-seq variants in that cell
  bam.path <- paste0("/home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/", ID, "/", DEPTH, "_recalibrated.bam")
  
  library(GenomicRanges)
  
  if(file.exists(bam.path) & file.exists(paths)){
    
    # 1. Read VCF
    vcf <- vcfR::read.vcfR(paths)
    
    # 2. Tidy & Format
    # Only select SNPs
    vcf <- vcfR::vcfR2tidy(vcf, format_fields = c("GT", "DP"))$fix %>%
      select(ChromKey,CHROM,POS,REF,ALT,QUAL,DP,AF1,VDB) %>%
      mutate(ID = ID,
             DEPTH = DEPTH,
             VARIANT = paste0(ChromKey, ":", POS, "-", REF, "/", ALT)) %>%
      filter(nchar(REF) == 1 & nchar(ALT) == 1) 
    
    # 3. Get Exome-seq variant expressed in this scRNA
    # From the mapping alignment in that cell, what is the exome-seq variant coverage ? 
    gr <- makeGRangesFromDataFrame(ref.pass, ignore.strand = TRUE, start.field = "POS", end.field = "POS", seqnames.field = "CHROM")
    bam.signal <- sapply(bamsignals::bamCoverage(bam.path, gr, verbose=FALSE), as.numeric)
    
    ref.pass$COVERAGE <- bam.signal
    ref.vars <- filter(ref.pass, COVERAGE > 2)$VARIANT
    
    # 4. Find Variant in reference exome-seq
    
    # 4.1. Filter out variants outside the exome-seq
    index <- suppressWarnings(findOverlaps(
      makeGRangesFromDataFrame(vcf, seqnames.field = "CHROM", start.field = "POS", end.field = "POS"),
      bed.exome)) %>% queryHits()
      
    # 4.2. Get variant
    # Currently only variants supported by >2 reads are considered
    vcf.dp <- vcf[index,] %>% 
      filter(DP > 2) %>%
      mutate(total.variant = n()) %>%
      group_by(DP, DEPTH) %>%
      summarise(ID = unique(ID),
                # Total number of variants
                total.variant = unique(total.variant),
                # Number of variant at that DP
                variant.dp.total = n(),
                # Total number of expressed variants
                available.variants = length(ref.vars),
                # scRNA variants in ref (true positives)
                variant.dp.inRef = sum(VARIANT %in% ref.vars),
                # scRNA variants not in ref (false positives)
                variant.dp.NotinRef = sum(!VARIANT %in% ref.vars),
                # scRNA variants not in ref and not among the dictionary of variants from other cells 
                variant.dp.NotinRef.NotInDict = sum( (!VARIANT %in% ref.vars) & (!VARIANT %in% common.vars)),
                # scRNA variants QUAL>20 in ref (true positives)
                variant.dp.inRef.QUAL = sum(VARIANT[QUAL > 20] %in% ref.vars),
                # scRNA variants QUAL>20 not in ref (false positive)
                variant.dp.NotinRef.QUAL = sum(!VARIANT[QUAL > 20] %in% ref.vars),
                # scRNA variants QUAL>20 not in ref and not among the dictionary of variants from other cells  (false positive)
                variant.dp.NotinRef.NotInDict.QUAL = sum( (!VARIANT[QUAL > 20] %in% ref.vars) & (!VARIANT[QUAL > 20] %in% common.vars)),
                overlap.refpass = 100*variant.dp.inRef / available.variants) %>%
      ungroup() 
    
    return(vcf.dp)
  } else {
    message("Can't find files")
  }
}

# 6. Run the scRNA VCF vs exome-seq comparison in parallel
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(doParallel))
threads <- 10

cl <- makeCluster(threads)
registerDoSNOW(cl)
clusterCall(cl, function() library(tidyverse))

vcf.dp <- foreach(i = 1:nrow(vcf.sp)) %dopar% {
  processVCF(paths = vcf.sp$PATHS[i], ID = vcf.sp$ID[i], DEPTH = vcf.sp$DEPTH[i], ref.pass = ref, bed.exome=bed.exome, common.vars = vars.dict$VARIANT[vars.dict$n > 10])
}

stopCluster(cl)

write_rds(vcf.dp, "/home/vincent.hahaut/data_storage/ORGANOIDS/vcf.results.rds")


# 7. Analyse the results
vcf.dp <- read_rds("/home/vincent.hahaut/data_storage/ORGANOIDS/vcf.results.rds")

# 7.0. Aggregate results
library(Seurat)
cellID_retained <- str_replace("^X", "", string = colnames(read_rds("/home/vincent.hahaut/data_storage/ORGANOIDS/ORGANOIDS.all.ssce.rds")))

vcf.dp.aggr <- vcf.dp %>% 
  bind_rows() %>%
  filter(ID %in% cellID_retained) %>%
  group_by(ID, DEPTH) %>%
  summarise(available.variants = unique(available.variants),
            total.vars = unique(total.variant),
            false_negative = unique(available.variants)-sum(variant.dp.inRef, na.rm = T),
            true_positive = sum(variant.dp.inRef, na.rm = T),
            true_positive.QUAL = sum(variant.dp.inRef.QUAL, na.rm = T),
            false_positive = sum(variant.dp.NotinRef, na.rm = T),
            false_positive.QUAL = sum(variant.dp.NotinRef.QUAL, na.rm = T),
            false_positive.filter = sum(variant.dp.NotinRef.NotInDict.QUAL, na.rm = T),
            false_positive.filter.QUAL = sum(variant.dp.NotinRef.NotInDict, na.rm = T),
            percentage_true = 100*sum(variant.dp.inRef, na.rm = T)/available.variants,
            percentage_true.QUAL = 100*sum(variant.dp.inRef.QUAL, na.rm = T)/available.variants
            ) %>%
  ungroup() %>%
  select(-ID) %>%
  pivot_longer(-DEPTH) %>%
  group_by(DEPTH, name) %>%
  summarise(m = mean(value),
            s = sd(value))

# 7.1. False positives, true positives, total number of variants, ...
cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

p.variants <- plot_grid(nrow = 2,
                        # available variants
                        ggplot(filter(vcf.dp.aggr, name %in% c("available.variants", "total.vars")), aes(x = as.numeric(DEPTH), fill = name, color = name)) + 
                          geom_hline(yintercept = seq(0,3500,500), linetype = "dashed", color = "darkgrey") +
                          geom_errorbar(aes(ymin = m - s, ymax = m + s), position = position_dodge(1000)) +
                          geom_line(aes(y = m, group = name), position = position_dodge(1000)) +
                          geom_point(aes(y = m), position = position_dodge(1000)) +
                          scale_color_manual(values = cols[4:5]) +
                          theme_cowplot(font_size = 16) + 
                          theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
                          scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                                             labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
                          scale_y_continuous(breaks = seq(0,4000,500), limits = c(0,4000)) +
                          xlab("Downsampled Reads") + ylab("SNPs") + ggtitle("Total SNPs"),
                        # False negative / True positive
    ggplot(filter(vcf.dp.aggr, name %in% c("true_positive.QUAL", "false_negative", "true_positive")), aes(x = as.numeric(DEPTH), fill = name, color = name)) + 
    geom_hline(yintercept = seq(0,1500,500), linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), position = position_dodge(1000)) +
    geom_line(aes(y = m, group = name), position = position_dodge(1000)) +
    geom_point(aes(y = m), position = position_dodge(1000)) +
    scale_color_manual(values = cols[1:3]) +
    theme_cowplot(font_size = 16) + 
    theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                       labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
    scale_y_continuous(breaks = seq(0,2000,500), limits = c(0,2000)) +
    xlab("Downsampled Reads") + ylab("SNPs") + ggtitle("True Variants"),
    # False Positive
  ggplot(filter(vcf.dp.aggr, name %in% c("false_positive", "false_positive.QUAL")), aes(x = as.numeric(DEPTH), fill = name, color = name)) + 
    geom_hline(yintercept = seq(0,2000,500), linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), position = position_dodge(1000)) +
    geom_line(aes(y = m, group = name), position = position_dodge(1000)) +
    geom_point(aes(y = m), position = position_dodge(1000)) +
    scale_color_manual(values = cols[6:8]) +
    theme_cowplot(font_size = 16) + 
    theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                       labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
    scale_y_continuous(breaks = seq(0,2500,500), limits = c(0,2500)) +
    xlab("Downsampled Reads") + ylab("SNPs") + ggtitle("False Positives"),
  # False Positive
  ggplot(filter(vcf.dp.aggr, name %in% c("percentage_true", "percentage_true.QUAL")), aes(x = as.numeric(DEPTH), fill = name, color = name)) + 
    geom_hline(yintercept = seq(0,60,20), linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), position = position_dodge(1000)) +
    geom_line(aes(y = m, group = name), position = position_dodge(1000)) +
    geom_point(aes(y = m), position = position_dodge(1000)) +
    scale_color_manual(values = cols[2:3]) +
    theme_cowplot(font_size = 16) + 
    theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                       labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
    scale_y_continuous(breaks = seq(0,80,20), limits = c(0,80)) +
    xlab("Downsampled Reads") + ylab("SNPs (%)") + ggtitle("True Positives (%)")
)

# 7.2. Influence of the variant coverage
vcf.dp.depth.aggr <- vcf.dp %>% 
  bind_rows() %>%
  mutate(bin = case_when(
    DP < 4 ~ "DP<3",
    DP == 4 ~ "DP=4",
    DP == 5 ~ "DP=5",
    DP > 5 & DP <= 10 ~ "5<DP≤10",
    DP > 10 & DP <= 20 ~ "11<DP≤20",
    DP > 20  ~ "DP>20"),
    bin = factor(bin, levels = c("DP<3","DP=4","DP=5","5<DP≤10","11<DP≤20","DP>20"))) %>%
  group_by(DEPTH, ID, bin) %>%
  summarise(true_positive = sum(variant.dp.inRef),
            false_positive = sum(variant.dp.NotinRef),
            true_positive.QUAL = sum(variant.dp.inRef.QUAL),
            false_positive.QUAL = sum(variant.dp.NotinRef.QUAL)) %>%
  ungroup() %>%
  select(-ID) %>%
  pivot_longer(cols = c(false_positive, true_positive, true_positive.QUAL, false_positive.QUAL)) %>%
  group_by(DEPTH, name, bin) %>%
  summarise(m = mean(value),
            s = sd(value))

p.depth <- plot_grid(nrow = 2,
  ggplot(filter(vcf.dp.depth.aggr, name %in% c("true_positive")), aes(x = as.numeric(DEPTH), fill = bin, color = bin)) + 
    geom_hline(yintercept = seq(0,500,100), linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), alpha = 0.5) +
    geom_line(aes(y = m, group = bin)) +
    geom_point(aes(y = m)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_cowplot(font_size = 20) + 
    theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
    scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                       labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
    scale_y_continuous(breaks = seq(0,400,100)) +
    xlab("Downsampled Reads") + ylab("Variants") + ggtitle("True Positives"),
  ggplot(filter(vcf.dp.depth.aggr, name %in% c("false_positive")), aes(x = as.numeric(DEPTH), fill = bin, color = bin)) + 
    geom_hline(yintercept = seq(0,700,100), linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = m - s, ymax = m + s), alpha = 0.5) +
    geom_line(aes(y = m, group = bin)) +
    geom_point(aes(y = m)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_cowplot(font_size = 20) +
    theme(legend.position = "none", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
    scale_x_continuous(breaks = c(5000, 10000,20000, 40000, 50000,75000,100000, 125000, 250000, 375000, 500000, 750000), 
                       labels = c("", "","", "", "50K","75K","100K", "125K", "250K", "375K", "500K", "750K")) +
    scale_y_continuous(breaks = seq(0,700,100)) +
    xlab("Downsampled Reads") + ylab("Variants") + ggtitle("False Positive")
)

ggsave(p.depth, bg = "white", filename = "/home/vincent.hahaut/data_storage/ORGANOIDS/FIGURES/ORGANOIDS_VARIANTS_DEPTH.tiff", dpi = 450, height = 8, width = 7)
ggsave(p.variants, bg = "white", filename = "/home/vincent.hahaut/data_storage/ORGANOIDS/FIGURES/ORGANOIDS_VARIANTS.tiff", dpi = 450, height = 8, width = 12)


# True positives / False positives relationship
tp_fp <- vcf.dp %>% 
  bind_rows() %>%
  group_by(ID, DEPTH) %>%
  summarise(available.variants = unique(available.variants),
            true_positive = sum(variant.dp.inRef),
            false_positive = sum(variant.dp.NotinRef)) %>%
  ungroup() %>%
  mutate(DEPTH = factor(DEPTH, levels = c("5000","10000","20000","40000","50000","75000","100000","125000","250000","375000","500000","750000")))

corrs <- ggplot() +
  geom_hline(yintercept = seq(1000,3000,1000), linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = seq(1000,3000,1000), linetype = "dashed", color = "darkgrey") +
  geom_point(data = tp_fp, aes(false_positive, true_positive, color = DEPTH), size = 0.75, alpha = 0.75) +
  geom_smooth(data = tp_fp, aes(false_positive, true_positive), method = lm, color = "black") +
  scale_color_manual(values = mycolors) +
  ylab("True Positives") + 
  xlab("False Positives") +
  theme_cowplot(font_size = 20) +
  geom_text(data = data.frame(x=2000,y=2800, DEPTH = "5000"), aes(x = 2000, y = 2800, label = "p-value < 2e-16"), size = 7) +
  geom_text(data = data.frame(x=2000,y=2800, DEPTH = "5000"), aes(x = 2000, y = 2600, label = "cor: 0.988599"), size = 7)


ggsave(corrs, bg = "white", filename = "/home/vincent.hahaut/data_storage/ORGANOIDS/FIGURES/ORGANOIDS_VARIANTS_corrs.tiff", dpi = 450, height = 6, width = 8)


####### 10x Variant Calling ############


# Resample 100 rods 10-times (replacement = F)
# Resample 10 rods 10-times (replacement = F)

# 1. Create the Sample Sheet
vcf.paths <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS_W18/VARIANTS/", pattern = ".vcf", recursive = TRUE, full.names = TRUE)
vcf.sp <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS_W18/VARIANTS/", pattern = ".vcf", recursive = TRUE, full.names = FALSE) 

pattern <- "GATK"
vcf.sp <- vcf.sp[grepl(vcf.paths, pattern = paste0("\\.", pattern, "\\."))]
vcf.paths <- vcf.paths[grepl(vcf.paths, pattern = paste0("\\.", pattern, "\\."))]

vcf.sp <- data_frame(PATHS = vcf.paths,
                     ID = str_split(vcf.sp, "/", simplify = TRUE)[,1],
                     ID_original = vcf.sp,
                     DEPTH = str_split(vcf.sp, "\\.", simplify = T)[,3])

ref.pass <- ref

# 2. Read the vcf
results <- data_frame()
for(i in 1:nrow(vcf.sp)){

  ID <- vcf.sp$ID[i]
  
  # 2. Get Bam info to filter out unexpressed exome-seq variants
  bam.path <- paste0("/home/vincent.hahaut/data_storage/ORGANOIDS_W18/VARIANTS/", ID, "/_recalibrated.bam")
  params <- Rsamtools::ScanBamParam(what=c("qname")) 
  bam <- Rsamtools::scanBam(bam.path, param = params)
  
  # 2.1. Read VCF
  vcf <- vcfR::read.vcfR(vcf.sp$PATHS[i])
  
  # 2.2. Tidy & Format
  # Only select SNPs
  vcf <- vcfR::vcfR2tidy(vcf, format_fields = c("GT", "DP"))$fix %>%
    dplyr::select(ChromKey,CHROM,POS,REF,ALT,QUAL,DP,AF1,VDB) %>%
    mutate(ID = ID,
           VARIANT = paste0(ChromKey, ":", POS, "-", REF, "/", ALT)) %>%
    filter(nchar(REF) == 1 & nchar(ALT) == 1) 
  
  # 2.3. Filter out variants in regions not covered by the exome-seq
  # 3. Get Exome-seq variant expressed in this scRNA
  # From the mapping alignment in that cell, what is the exome-seq variant coverage ? 
  gr <- GenomicRanges::makeGRangesFromDataFrame(ref.pass, ignore.strand = TRUE, start.field = "POS", end.field = "POS", seqnames.field = "CHROM")
  bam.signal <- sapply(bamsignals::bamCoverage(bam.path, gr, verbose=FALSE), as.numeric)
  
  ref.pass$COVERAGE <- bam.signal
  ref.vars <- filter(ref.pass, COVERAGE > 2)$VARIANT
  
  # 4. Find Variant in reference exome-seq
  
  # 4.1. Filter out variants outside the exome-seq
  index <- suppressWarnings(findOverlaps(
    makeGRangesFromDataFrame(vcf, seqnames.field = "CHROM", start.field = "POS", end.field = "POS"),
    bed.exome)) %>% queryHits()
  
  # 5. Overlap 10x
  results <- bind_rows(results, 
    vcf[index,] %>% 
    filter(DP > 2) %>%
    mutate(total.variant = n()) %>%
    summarise(ID = unique(ID),
              # Total number of variants
              total.variant = unique(total.variant),
              # Total number of expressed variants
              available.variants = length(ref.vars),
              # scRNA variants in ref (true positives)
              variant.dp.inRef = sum(VARIANT %in% ref.vars),
              # scRNA variants not in ref (false positives)
              variant.dp.NotinRef = sum(!VARIANT %in% ref.vars),
              overlap.refpass = 100*variant.dp.inRef / available.variants,
              # Extra info
              nReads = length(bam[[1]]$qname),
              nCells = ifelse(grepl("_10_", ID), 10, 100))
  )
}


# 3. Aggregate results
results <- results %>%
  mutate(true_positive = variant.dp.inRef,
         false_positive = variant.dp.NotinRef,
         false_negative = available.variants-variant.dp.inRef)

# 4. Visualize results
cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

p.variants.10x <- plot_grid(
  results %>%
    dplyr::select(nCells, true_positive, false_positive, false_negative) %>%
    pivot_longer(-nCells) %>%
    mutate(name = factor(name, levels = c("true_positive", "false_positive", "false_negative")),
           nCells = paste0("nCells:", nCells)) %>%
    ggplot(aes(x = name, y = value, fill = name)) + 
    geom_hline(yintercept = seq(0,1000,500), linetype = "dashed", color = "darkgrey") +
    facet_wrap("nCells") + 
    geom_violin() +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    scale_fill_manual(values = cols[c(2,6,1)]) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(breaks = seq(0,1000,250)) +
    xlab("") + ylab("Variants"),
  results %>%
    dplyr::select(nCells, total.variant, available.variants) %>%
    pivot_longer(-nCells) %>%
    mutate(name = factor(name, levels = c("available.variants", "total.variant")),
           nCells = paste0("nCells:", nCells)) %>%
    ggplot(aes(x = name, y = value, fill = name)) + 
    geom_hline(yintercept = seq(0,1000,500), linetype = "dashed", color = "darkgrey") +
    facet_wrap("nCells") + 
    geom_violin() +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    scale_fill_manual(values = cols[c(4,5)]) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(breaks = seq(0,1000,250)) +
    xlab("") + ylab("Variants")
)

results %>%
  group_by(nCells) %>%
  summarise(mean(nReads), sd(nReads))

ggsave(p.variants.10x, filename = "/home/vincent.hahaut/data_storage/ORGANOIDS/FIGURES/10x_variant_calling.tiff", dpi = 450, height = 6, width = 10, bg = "white")

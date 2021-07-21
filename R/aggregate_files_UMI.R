# All samples involved in the FS-UMI comparisons where processed with 4_SINGLE_END_UMI.sh and regrouped in:
path <- "/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/"
out.prefix <- 

# 1. Internal / 5' UMI gene detection
if(!file.exists("/home/vincent.hahaut/data_storage/FS_intermediate_file/UMI.txt")){
  
  ID <- list.files(path)
  
  df.umi <- df.all <- data_frame()
  for(i in 1:length(ID)){
    path.id <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID[i], "/DOWNSAMPLE/")
    
    # Downsamplings
    for(down in c("10000","25000","50000","75000","100000", "125000", "250000")){
      path.id.down <- paste0(path.id, "DOWN_", down)
      
      # UMI or internal or both reads
      for(type in c("UMI", "ALL", "INTERNAL")){
        
        if(type %in% c("ALL", "INTERNAL")){
          path.id.down.exon <- paste0(path.id.down, "/FEATURECOUNTS/", ID[i], "_", type, "_", down, ".feature.EXON.txt")
          path.id.down.gene <- paste0(path.id.down, "/FEATURECOUNTS/", ID[i], "_", type, "_", down, ".feature.GENE.txt")
          
          if(file.exists(path.id.down.exon)){
            exon <- as.vector(unlist(read_tsv(path.id.down.exon, comment = "#", col_names = TRUE)[,7]))
            gene <- as.vector(unlist(read_tsv(path.id.down.gene, comment = "#", col_names = TRUE)[,7]))
            
            df.all <- bind_rows(df.all,
                                data_frame(ID = ID[i],
                                           TYPE = type,
                                           DOWN = down,
                                           count.0.exon = sum(exon > 0),
                                           count.1.exon = sum(exon > 1),
                                           count.5.exon = sum(exon > 5),
                                           count.0.gene = sum(gene > 0),
                                           count.1.gene = sum(gene > 1),
                                           count.5.gene = sum(gene > 5))
            )
            
          }
          
        } else if(type == "UMI"){
          path.id.down.umi <- paste0(path.id.down, "/FEATURECOUNTS/umi.counts.", down, ".tsv.gz")
          if(file.exists(path.id.down.umi)){
            umi <- as.vector(unlist(read_tsv(path.id.down.umi, comment = "#", col_names = TRUE)[,2]))
            
            df.umi <- bind_rows(df.umi,
                                data_frame(ID = ID[i],
                                           TYPE = type,
                                           DOWN = down,
                                           MOLECULES = sum(umi),
                                           count.0.umi = sum(umi > 0),
                                           count.1.umi = sum(umi > 1),
                                           count.5.umi = sum(umi > 5)
                                )
            )
            
          }
          
        }
      }
    }
  }
  write.table(df.umi, "/home/vincent.hahaut/data_storage/FS_intermediate_file/UMI.txt", sep = "\t", row.names = F, quote = F,)
  write.table(df.all, "/home/vincent.hahaut/data_storage/FS_intermediate_file/ALL.txt", sep = "\t", row.names = F, quote = F,)
}

####################################################
# 2. Proportion of UMI reads vs ALL (UMI+INTERNAL) #
####################################################

if(!file.exists("/home/vincent.hahaut/data_storage/FS_intermediate_file/reads.umi.txt")){
  
  ID <- list.files(path)
  path.all <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/FASTQ/allreads.R1.fq.gz")
  path.interne <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/FASTQ/internal.R1.fq.gz")
  path.umi <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/FASTQ/umi.R1.fq.gz")
  
  reads.all <- lapply(1:length(path.all), function(x) data_frame(ALL = ShortRead::readFastq(path.all[x]) %>% length, ID = ID[x]))
  reads.umi <- lapply(1:length(path.umi), function(x) data_frame(ALL = ShortRead::readFastq(path.umi[x]) %>% length, ID = ID[x]))
  reads.all <- bind_rows(reads.all) %>% mutate(TYPE = "ALL")
  reads.umi <- bind_rows(reads.umi) %>% mutate(TYPE = "UMI")
  
  write.table(reads.all, "/home/vincent.hahaut/data_storage/FS_intermediate_file/reads.all.txt", sep = "\t", row.names = F, quote = F)
  write.table(reads.umi, "/home/vincent.hahaut/data_storage/FS_intermediate_file/reads.umi.txt", sep = "\t", row.names = F, quote = F)
  
}

#########################
# 3. Gene Body Coverage #
#########################

ID <- list.files(path)
path.interne <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/ReSQC/", ID, "_INTERNAL_geneBody.all.geneBodyCoverage.txt")
path.umi <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/ReSQC/", ID, "_UMI_geneBody.all.geneBodyCoverage.txt")

# see aggregate_files.R for getRESeQC_coverage
# 3.1. Get Coverage internal reads
df.internal <- getRESeQC_coverage(path.interne, ID, thread = 5) %>%
  mutate(TYPE = "INTERNAL",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3 Hagemann-J.",
           grepl("AAGCA", ID) ~ "TSO-AAGCA + STRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA + STRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC + STRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC + STRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC + STRT-dT",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC + STRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA + STRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI + STRT-dT",
           grepl("276_LA", ID) ~ "LA-TSO-CAGCA + STRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC + FS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA + FS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI + FS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA + SS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI + SS3-dT")) %>%
  filter(GROUP %in% condition)

# 3.2. Get coverage UMI reads
df.umi <- getRESeQC_coverage(path.umi, ID, thread = 5) %>%
  mutate(TYPE = "UMI",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3 Hagemann-J.",
           grepl("AAGCA", ID) ~ "TSO-AAGCA + STRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA + STRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC + STRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC + STRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC + STRT-dT",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC + STRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA + STRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI + STRT-dT",
           grepl("276_LA", ID) ~ "LA-TSO-CAGCA + STRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC + FS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA + FS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI + FS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA + SS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI + SS3-dT")) %>%
  filter(GROUP %in% condition)

# 3.3. Combine UMI / Internal reads
df.coverage <- bind_rows(df.internal, df.umi) %>%
  mutate(GROUP = 
           factor(GROUP, 
                  levels = c("SS3", "SS3 Hagemann-J.","TSO-ATAAC + STRT-dT","TSO-CAGCA + STRT-dT","TSO-FS-UMI + STRT-dT","TSO-ATAAC + FS-dT","TSO-CAGCA + FS-dT","TSO-FS-UMI + FS-dT","TSO-CAGCA + SS3-dT","TSO-FS-UMI + SS3-dT", "TSO-AAGCA + STRT-dT","TSO-CATCA + STRT-dT","TSO-CTGAC + STRT-dT","TSO-ATGAC + STRT-dT",
                             "TSO-CTAAC + STRT-dT", "LA-TSO-CAGCA + STRT-dT"))) 

########################
# 4. Read Distribution #
########################

ID <- list.files(path)
path.interne <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/ReSQC/", ID, "_INTERNAL_readDistribution.txt")
path.umi <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/ReSQC/", ID, "_UMI_readDistribution.txt")

# 4.1. Internal Reads
df.internal <- getRESeQC_readDistribution(path.interne, ID, thread = 5) %>%
  mutate(TYPE = "INTERNAL",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3\nHagemann-J.",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC\nSTRT-dT",
           grepl("AAGCA", ID) ~ "TSO-AAGCA\nSTRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA\nSTRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC\nSTRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC\nSTRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC\nSTRT-dT",
           grepl("276_LA", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI\nSTRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC\nFS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA\nFS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI\nFS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA\nSS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI\nSS3-dT")) %>%
  filter(GROUP %in% condition)

# 4.2. UMI Reads
df.umi <- getRESeQC_readDistribution(path.umi, ID, thread = 5) %>%
  mutate(TYPE = "UMI",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3\nHagemann-J.",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC\nSTRT-dT",
           grepl("AAGCA", ID) ~ "TSO-AAGCA\nSTRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA\nSTRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC\nSTRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC\nSTRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC\nSTRT-dT",
           grepl("276_LA", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI\nSTRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC\nFS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA\nFS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI\nFS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA\nSS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI\nSS3-dT")) %>%
  filter(GROUP %in% condition)

# 4.3. Combine UMI / Internal
df.readDistrib <- bind_rows(df.internal, df.umi) %>%
  dplyr::select(-ID) %>%
  gather(distribution, val, -GROUP, -TYPE) %>%
  filter(!distribution %in% c("TSS_up_1kb", "TSS_up_5kb", "TSS_down_1kb", "TSS_down_1kb",  "TES_down_10kb", "TES_down_1kb" , "TES_down_5kb", "TSS_up_10kb", "totalReadTags")) %>%
  mutate(distribution = case_when(distribution == "UTR5_Exons" ~ "5'UTR Exons",
                                  distribution == "UTR3_Exons" ~ "3'UTR Exons",
                                  distribution == "CDS_Exons" ~ "CDS Exons",
                                  distribution == "Introns" ~ "Introns",
                                  distribution == "Intergenic" ~ "Intergenic"),
         distribution = factor(distribution, levels = rev(c("CDS Exons", "Introns", "5'UTR Exons", "3'UTR Exons","Intergenic")))) %>%
  # Compute the per-group % reads
  group_by(GROUP, TYPE) %>%
  mutate(total = sum(val)) %>%
  ungroup() %>%
  group_by(GROUP, TYPE, distribution) %>%
  summarise(perc = 100* sum(val) / unique(total)) %>%
  ungroup() %>%
  mutate(facet = case_when(grepl("SS3$|^SS3", GROUP) ~ "SS3",
                           grepl("TSO-FS-UMI", GROUP) ~ "FS-UMI",
                           TRUE  ~ "SPACER"))
  


#########################
# 5. STAR Mapping Stats #
#########################


ID <- list.files(path)
path.interne <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/STAR/", ID, "_INTERNAL_Log.final.out")
path.umi <- paste0("/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/", ID, "/STAR/", ID, "_UMI_Log.final.out")

# 5.1. Read the files
df.internal <- getSTARLog(path.interne, ID, thread = 5) %>%
  mutate(TYPE = "INTERNAL",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3\nHagemann-J.",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC\nSTRT-dT",
           grepl("AAGCA", ID) ~ "TSO-AAGCA\nSTRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA\nSTRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC\nSTRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC\nSTRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC\nSTRT-dT",
           grepl("276_LA", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI\nSTRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC\nFS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA\nFS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI\nFS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA\nSS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI\nSS3-dT")) %>%
  filter(GROUP %in% condition)

df.umi <- getSTARLog(path.umi, ID, thread = 5) %>%
  mutate(TYPE = "UMI",
         GROUP = case_when(
           grepl("5_SS3|6_SS3", ID) ~ "SS3",
           grepl("SS3fwd", ID) ~ "SS3\nHagemann-J.",
           grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC\nSTRT-dT",
           grepl("AAGCA", ID) ~ "TSO-AAGCA\nSTRT-dT",
           grepl("CATCA", ID) ~ "TSO-CATCA\nSTRT-dT",
           grepl("CTGAC", ID) ~ "TSO-CTGAC\nSTRT-dT",
           grepl("ATGAC", ID) ~ "TSO-ATGAC\nSTRT-dT",
           grepl("CTAAC", ID) ~ "TSO-CTAAC\nSTRT-dT",
           grepl("276_LA", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA\nSTRT-dT",
           grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI\nSTRT-dT",
           grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC\nFS-dT",
           grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA\nFS-dT",
           grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI\nFS-dT",
           grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA\nSS3-dT",
           grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI\nSS3-dT")) %>%
  filter(GROUP %in% condition)

# 5.2. Aggregate the results
df.star <- bind_rows(df.internal, df.umi) %>%
  dplyr::select(ID, N_raw_reads, N_uniquely_mapped_reads, P_too_short_read, P_multiple_Loci, GROUP, TYPE) %>% 
  mutate(N_raw_reads = as.numeric(N_raw_reads), 
         N_uniquely_mapped_reads = as.numeric(N_uniquely_mapped_reads),  
         P_too_short_read = as.numeric(P_too_short_read), 
         P_multiple_Loci = as.numeric(P_multiple_Loci)) %>%
  group_by(ID) %>%
  mutate(N_too_short_read = P_too_short_read/100 * N_raw_reads,
         N_multiple_Loci = P_multiple_Loci/100 * N_raw_reads) %>%
  dplyr::select(-P_too_short_read, -N_raw_reads, -P_multiple_Loci) %>%
  gather(reads, n, -ID, -GROUP, -TYPE) %>%
  mutate(reads = case_when(
    reads == "N_uniquely_mapped_reads" ~ "Uniquely Mapped",
    reads == "N_too_short_read" ~ "Unmapped",
    reads == "N_multiple_Loci" ~ "Multiple Loci"),
    reads = factor(reads, levels = rev(c("Uniquely Mapped", "Multiple Loci", "Unmapped")))) %>%
  # Compute the per-group % for barplots
  group_by(GROUP, TYPE) %>%
  mutate(total = sum(n)) %>%
  group_by(GROUP, TYPE, reads) %>%
  summarise(mp = 100*sum(n)/unique(total)) %>%
  ungroup() %>%
  mutate(facet = case_when(grepl("SS3$|^SS3", GROUP) ~ "SS3",
                           grepl("TSO-FS-UMI", GROUP) ~ "FS-UMI",
                           TRUE  ~ "SPACER")) 



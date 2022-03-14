#' createGTFfile
#' 
#' Create a GTF containing the reduced exon-intron positions
#' Adapted from Devon Ryan (https://www.biostars.org/p/165226/) - parallelised and allow genes with no introns.
#' @param gtf.path, GTF file from gencode or ensembl
#' @param prefix, output path prefix
#' @param threads, number of threads to use for the intron search

createGTFfile <- function(gtf.path = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/...gtf", prefix = NULL, threads = 30){
  
  # 0. Libraries
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(doSNOW))
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(tidyverse))
  
  # 1. Read the GTF file
  print("1. Read GTF")
  gtf <- makeTxDbFromGFF(gtf.path)
  
  # 2. Get the exonic sequences and reduce the overlap per gene.
  print("2. Get the exons")
  exons <- exonsBy(gtf, by="gene")
  exons <- IRanges::reduce(exons)
  
  # 3. Get the intronic sequences
  print("3. Get the introns")
  
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  clusterCall(cl, function() library(GenomicFeatures))
  
  res <- foreach(i = 1:length(exons)) %dopar% {
    x <- exons[[i]]
    if(length(x) > 1){
      gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
                                                           end=max(end(x))), 
                   strand=strand(x)[1])
      db = disjoin(c(x, gr))
      ints = db[countOverlaps(db, x) == 0]
      ints$ID <- names(exons)[i]
      #Add an ID
      if(as.character(strand(ints)[1]) == "-") {
        ints$exon_id = c(length(ints):1)
      } else {
        ints$exon_id = c(1:length(ints))
      }
      return(ints)
    }
  }
  stopCluster(cl)
  
  # 4. Tidy and regroup the results
  print("4. Regroup")
  introns <- lapply(res, as.data.frame) %>%
    bind_rows() %>%
    mutate(INTRON_ID = paste0(ID, "_", exon_id)) %>% 
    ungroup() %>%
    mutate(TYPE = "INTRON")
  exons <- as.data.frame(exons) %>% 
    group_by(group_name) %>% 
    mutate(exon_num = 1:n()) %>% 
    dplyr::select(-group) %>% 
    dplyr::rename("ID" = group_name) %>% 
    mutate(EXON_ID = paste0(ID, "_",exon_num)) %>% 
    ungroup() %>%
    mutate(TYPE = "EXON")
  
  introns <- dplyr::rename(introns, intron_num = "exon_id")
  
  # 5. Add Gene Names
  gtf.open <- rtracklayer::import(gtf.path) %>% as_tibble()
  introns <- left_join(introns, gtf.open %>% dplyr::select(gene_id, gene_name) %>% distinct(), by = c("ID" = "gene_id"))
  exons <- left_join(exons, gtf.open %>% dplyr::select(gene_id, gene_name) %>% distinct(), by = c("ID" = "gene_id"))
  
  # 6. Create a GTF file
  
  # 6.1. Introns
  print("5. To GTF")
  intron.gtf <- data.frame(
    seqnames = introns$seqnames,
    source = ".",
    feature = "intron",
    start = introns$start,
    end = introns$end,
    score = ".",
    strand = introns$strand,
    frame = ".",
    attribute = paste0('gene_id "', introns$ID, '"; ', 
                       'gene_name "', introns$gene_name, '"; ', 
                       'intron_num "', introns$intron_num, '"; ', 
                       'intron_id "', introns$INTRON_ID, '"')) 
  
  # 6.2. Exons
  exon.gtf <- data.frame(
    seqnames = exons$seqnames,
    source = ".",
    feature = "exon",
    start = exons$start,
    end = exons$end,
    score = ".",
    strand = exons$strand,
    frame = ".",
    attribute = paste0('gene_id "', exons$ID, '"; ', 
                       'gene_name "', exons$gene_name, '"; ', 
                       'exon_num "', exons$exon_num, '"; ', 
                       'exon_id "', exons$EXON_ID, '"')) 
  
  # 6.3. Unified reference
  introns$featureID <- paste0("INTRON-", introns$INTRON_ID, "_", introns$strand)
  exons$featureID <- paste0("EXON-", exons$EXON_ID, "_", exons$strand)
  introns$feature_name <- paste0("INTRON-", introns$gene_name, "_", introns$strand)
  exons$feature_name <- paste0("EXON-", exons$gene_name, "_", exons$strand)
  introns$gene_id <- str_replace(introns$INTRON_ID, "_.+$", "")
  exons$gene_id <- str_replace(exons$EXON_ID, "_.+$", "")
    
  introns$INTRON_ID <- NULL
  introns$intron_num <- NULL
  exons$EXON_ID <- NULL
  exons$exon_num <- NULL
  
  unified <- bind_rows(exons, introns)
  unified <- data.frame(seqnames = unified$seqnames,
                        source = ".",
                        feature = "all",
                        start= unified$start,
                        end = unified$end,
                        score = ".",
                        strand = unified$strand,
                        frame = ".",
                        attribute = paste(paste0("feature=",unified$featureID),
                                          paste0("feature_name=", unified$feature_name),
                                          paste0("gene_id=", unified$gene_id),
                                          paste0("gene_name=", unified$gene_name), sep = "; ")) %>%
    arrange(seqnames, start, end) %>%
    as.data.frame()
  
  
  # 7. Save as GTF
  print("6. Save data")
  write.table(intron.gtf, paste0(prefix, ".introns.gtf"), sep = "\t", col.names = F, row.names = F, quote = F)
  write.table(exon.gtf, paste0(prefix, ".exons.gtf"), sep = "\t", col.names = F, row.names = F, quote = F)
  write.table(unified, paste0(prefix, ".exonsANDintrons.gtf"), sep = "\t", col.names = F, row.names = F, quote = F)
  
}

createGTFfile(gtf.path = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf", prefix = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly", threads = 30)

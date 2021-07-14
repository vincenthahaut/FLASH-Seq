#' checkInputPath
#' 
#' Check if provided path exists and remove empty links
#' @param paths = vector of path
#' @param sample.ids = vector of sample ids
checkInputPath <- function(paths = NULL, sample.ids = NULL){
  
  exists <- file.exists(paths)
  message(paste0("The Following Files Do Not Exists: ", paths[!exists]))
  
  paths <- paths[exists]
  sample.ids <- sample.ids[exists]
  
  return(list(paths = paths, sample.ids = sample.ids))
}

######################
# Mapping Statistics # 
######################

#' getSTARLog
#' 
#' Get STAR mapping statistics from the log.final.
#' @param sample.ids = vector of sample.ids.
#' @param paths = vector of path/to/STAR_log.final.out
#' @param threads = Number of threads for parallel processing.
getSTARLog <- function(paths = NULL, sample.ids = NULL, threads = 5){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(doSNOW))
  suppressPackageStartupMessages(library(doParallel))
  if(threads >= detectCores()){stop("Requested too many cores")}
  
  # 1. Check if provided paths exist
  trueFiles <- checkInputPath(paths, sample.ids)
  paths <- trueFiles$paths
  sample.ids <- trueFiles$sample.ids
  
  # 2. Process All Files in Parallel
  message("1. Prepare Multi-threading")
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  clusterCall(cl, function() library(tidyverse))
  
  pb <- txtProgressBar(max = length(paths), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  message("2. Start Reading")
  
  STAR_stats <- foreach(i = 1:length(paths), .options.snow = opts) %dopar% {
    
    read.delim(paths[i], skip = 3, sep = "\t", header = FALSE) %>%
      dplyr::select(V2) %>%
      mutate(V2 = trimws(str_replace_all(V2, "%", ""))) %>%
      dplyr::slice(c(2,3,5,6,7,8,11,20,21,22,23,25,27,28,29)) %>%
      mutate(info = c("N_raw_reads", 
                      "Average_read_length", 
                      "N_uniquely_mapped_reads",  
                      "P_uniquely_mapped_reads",
                      "Average_mapped_read_length", 
                      "N_splice", 
                      "N_splice_GCAG",
                      "N_multiple_Loci", 
                      "P_multiple_Loci", 
                      "N_many_Loci", 
                      "P_many_Loci", 
                      "N_unmapped_mismatch",
                      "N_too_short_read", 
                      "P_too_short_read", 
                      "N_unmapped_other")) %>%
      spread(info, V2) %>%
      mutate(N_umapped_all = sum(as.numeric(N_too_short_read), as.numeric(N_unmapped_other), as.numeric(N_unmapped_mismatch)),
             N_multi = sum(as.numeric(N_multiple_Loci), as.numeric(N_many_Loci))) %>%
      dplyr::select(-N_too_short_read, -N_unmapped_other, -N_unmapped_mismatch, -N_multiple_Loci, -N_many_Loci) %>%
      mutate(ID = sample.ids[i])
  }
  
  close(pb)
  stopCluster(cl)
  
  # 3. Regroup
  message("3. Aggregate Results")
  bind_rows(STAR_stats) %>%
    dplyr::select(ID, everything()) %>%
    return()
  
}

#' getRESeQC_coverage
#' 
#' Get the gene coverage as measured by ReSQC.
#' @param sample.ids = vector of sample.ids.
#' @param paths = vector of path/to/ReSQC_coverage.txt
#' @param threads = Number of threads for parallel processing.
getRESeQC_coverage <- function(paths = NULL, sample.ids = NULL, threads = 10){
  
  if(!require(data.table)){install.packages("data.table")}
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(doSNOW))
  suppressPackageStartupMessages(library(doParallel))
  if(threads >= detectCores()){stop("Requested too many cores")}
  
  # 1. Check if provided paths exist
  trueFiles <- checkInputPath(paths, sample.ids)
  paths <- trueFiles$paths
  sample.ids <- trueFiles$sample.ids
  
  # 2. Process All Files in Parallel
  message("1. Prepare Multi-threading")
  cl <- makeCluster(threads)
  registerDoSNOW(cl)
  clusterCall(cl, function() library(tidyverse))
  clusterCall(cl, function() library(data.table))
  
  message("2. Start Reading")
  
  pb <- txtProgressBar(max = length(paths), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  geneBodyCoverage <- foreach(i = 1:length(paths), .options.snow = opts) %dopar% {
    # Added if size > 0.
    # Testing gene length vs coverage some file were created empty. No need to pick them for now.
    # Remove later
    if(file.info(paths[i])$size > 0){
      data.table::fread(paths[i], header = TRUE, check.names = FALSE)[1,-1] %>%
        tidyr::as_tibble() %>%
        dplyr::mutate(ID = sample.ids[i])
    }
  }
  
  close(pb)
  stopCluster(cl)
  
  # 3. Regroup
  message("3. Aggregate Results")
  bind_rows(geneBodyCoverage) %>%
    dplyr::select(ID, everything()) %>%
    return()
  
}

#' getRESeQC_readDistribution
#' 
#' Get the intron / exon / intergenic distribution of the reads in a BAM files (according to ReSQC).
#' @param sample.ids = vector of sample.ids.
#' @param paths = vector of path/to/ReSQC_readDistribution.txt
#' @param threads = Number of threads for parallel processing.
getRESeQC_readDistribution <- function(paths = NULL, sample.ids = NULL, threads = NULL){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(doSNOW))
  suppressPackageStartupMessages(library(doParallel))
  if(threads >= detectCores()){stop("Requested too many cores")}
  
  # 1. Check if provided paths exist
  trueFiles <- checkInputPath(paths, sample.ids)
  paths <- trueFiles$paths
  sample.ids <- trueFiles$sample.ids
  
  # 2. Process All Files in Parallel
  message("1. Prepare Multi-threading")
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  clusterCall(cl, function() library(tidyverse))
  
  message("2. Start Reading")
  pb <- txtProgressBar(max = length(paths), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  readDistrib <- foreach(i = 1:length(paths), .options.snow = opts) %dopar% {
    if(file.info(paths[i])$size > 0){
      data.table::fread(paths[i], header = TRUE, sep = " ", skip = 3) %>%
        select(Group, Tag_count) %>%
        mutate(ID = sample.ids[i]) %>%
        spread(Group, Tag_count) %>%
        dplyr::rename("UTR3_Exons" = `3'UTR_Exons`,
                      "UTR5_Exons" = `5'UTR_Exons`) %>%
        mutate(totalReadTags = as.numeric(data.table::fread(paths[i], header = FALSE, sep = " ", fill=TRUE, select = 1:3)$V3[2])) %>%
        rowwise() %>%
        mutate(Intergenic = totalReadTags-UTR3_Exons-UTR5_Exons-CDS_Exons-Introns-TSS_up_10kb,TES_down_10kb)
    } 
    
  }
  
  close(pb)
  stopCluster(cl)
  
  # 3. Regroup
  message("3. Aggregate Results")
  bind_rows(readDistrib) %>%
    dplyr::select(ID, everything()) %>%
    return()
}

###############
# Gene Counts # 
###############


#' getFeatureCounts
#' 
#' Read FeatureCount files in parallel. 
#' Assumes that all the files have the same dimension (= same GTF).
#' @param paths = vector of sample paths
#' @param sample.ids = vector of sample ids
#' @param mode = NULL or "downsampling". If downsampling is chosen, will attempt to compute the median number of genes detected per downsampling for each sample. Expects "ID_1", "ID_2", "ID_3", ... naming scheme to find the downsampled replicates. 
#' @param threads = Number of threads for parallel processing (Only for reading).
#' @param rowRange, boolean. Should you or not add rowRanges info ?
getFeatureCount <- function(paths = NULL, sample.ids = NULL, threads = 20, mode = NULL, rowRange = TRUE){
  
  if(!require(Rsamtools)){BiocManager::install("SingleCellExperiment")}
  if(!require(data.table)){library("data.table")}
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(doParallel))
  if(threads >= detectCores()){stop("Requested too many cores")}
  
  # 1. Check if provided paths exist
  trueFiles <- checkInputPath(paths, sample.ids)
  paths <- trueFiles$paths
  sample.ids <- trueFiles$sample.ids
  
  # 2. Read All Files in Parallel
  message("1. Prepare for multi-threading")
  geneID <- data.table::fread(paths[1], header = TRUE, select = 1, skip = 1, sep = "\t")
  
  message("2. Read Files")
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  mycounts.all <- foreach(i = 1:length(paths)) %dopar% {
    data.table::fread(paths[i], header = TRUE, check.names = FALSE, select = 7, skip = 1, sep = "\t")
  }
  stopCluster(cl)
  
  # 3. Regroup
  message("3. Regroup data")
  mycounts.all <- bind_cols(mycounts.all)
  mycounts.all <- bind_cols(geneID, mycounts.all)
  colnames(mycounts.all) <- c("gene_id", sample.ids)
  
  # 4. Average counts in downsamplings
  if(mode == "downsampling"){
    message("Downsampling Average Requested - Slow Process")
    
    mycounts.tidy <- mycounts.all %>%
      gather(ID, count, -gene_id) 
    
    # NB: Assume that the last part of the sample id is "_1" to "_100"
    # Will get the median count for each
    mycounts.all <- mycounts.tidy %>% head(1e6) %>%
      mutate(replicate = str_replace(ID, paste0(paste0("_", 1:100, "$"), collapse = "|"), "")) %>%
      group_by(gene_id, replicate) %>%
      summarise(count.med = median(count)) %>%
      ungroup() %>%
      spread(replicate, count.med)
  }
  
  # 5. Gene Info
  # Collapse unique chr: Some genes are located on both chrY and chrX
  if(rowRange == TRUE){ 
    message("4. Add rowRanges information")
    geneInfo <- read.table(paths[1], header = TRUE, skip = 1)[,1:6]
    geneInfo <- geneInfo %>%
      summarise(seqname = sapply(str_split(geneInfo$Chr, ";"), function(x) paste0(unique(x[[1]]), collapse = ";")),
                start = sapply(str_split(geneInfo$Start, ";"), function(x) min(x[[1]])),
                end = sapply(str_split(geneInfo$End, ";"), function(x) max(x[[1]])),
                strand = sapply(str_split(geneInfo$Strand, ";"), function(x) paste0(unique(x[[1]]), collapse = ";")),
                gene_name = Geneid) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # 6. Returns a SingleCellExperiment object
    message("5. To SingleCellExperiment")
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.frame(mycounts.all[,-1], row.names = mycounts.all$gene_id)),
                                                      rowRanges = geneInfo)
  } else {
    message("5. To SingleCellExperiment")
    
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.frame(mycounts.all[,-1], row.names = geneID$Geneid)))
  }
  return(sce)
  
}


#' getHTSeq
#' 
#' Get the gene count from HTSeq data.
#' Assumes that all the files have the same dimension (= same GTF).
#' @param sample.id = ID of the sample
#' @param path = path/to/STAR/htseq.txt
#' @param mode = NULL or "downsampling". If downsampling is chosen, will attempt to compute the median number of genes detected per downsampling for each sample. Expects "ID_1", "ID_2", "ID_3", ... naming scheme to find the downsampled replicates. 
#' @param path = path/to/file.gtf. The GTF used to create these files. 
#' @param threads = Number of threads for parallel processing (Only for reading).
getHTSeq <- function(path = NULL, sample.id = NULL, gtf = NULL, threads = 20){
  
  
  if(!require(Rsamtools)){BiocManager::install("SingleCellExperiment")}
  if(!require(data.table)){library("data.table")}
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(doParallel))
  if(threads >= detectCores()){stop("Requested too many cores")}
  
  # 1. Check if provided paths exist
  trueFiles <- checkInputPath(paths, sample.ids)
  paths <- trueFiles$paths
  sample.ids <- trueFiles$sample.ids
  
  # 2. Read All Files in Parallel
  geneID <- data.table::fread(paths[1], header = FALSE, select = 1, sep = "\t")
  
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  mycounts.all <- foreach(i = 1:length(paths)) %dopar% {
    data.table::fread(paths[i], header = FALSE, select = 2, sep = "\t")
  }
  stopCluster(cl)
  
  # 3. Regroup
  mycounts.all <- bind_cols(mycounts.all)
  mycounts.all <- bind_cols(geneID, mycounts.all)
  colnames(mycounts.all) <- c("gene_id", sample.ids)
  
  mycounts.all <- mycounts.all[1:(grep("__no_feature", mycounts.all$gene_id)-1),]
  
  
  # 4. Average counts in downsamplings
  if(mode == "downsampling"){
    message("Downsampling Average Requested - Slow Process")
    
    mycounts.tidy <- mycounts.all %>%
      gather(ID, count, -gene_id) 
    
    # NB: Assume that the last part of the sample id is "_1" to "_100"
    # Will get the median count for each
    mycounts.all <- mycounts.tidy %>% head(1e6) %>%
      mutate(replicate = str_replace(ID, paste0(paste0("_", 1:100, "$"), collapse = "|"), "")) %>%
      group_by(gene_id, replicate) %>%
      summarise(count.med = median(count)) %>%
      ungroup() %>%
      spread(replicate, count.med)
  }
  
  # 5. Gene Info
  # Collapse unique chr: Some genes are located on both chrY and chrX
  gtf <- rtracklayer::import(gtf)
  
  # 6. Returns a SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.frame(mycounts.all[,-1], row.names = mycounts.all$gene_id)),
                                                    rowRanges = gtf)
  
  return(sce)
  
}

#' getUMICounts
#' 
#' Extract count information from umi_tools count files. 
#' @param sample.ids vector of ids.
#' @param paths vector of /path/to/featurecounts/umi.counts.tsv.gz
#' @return SingleCellExperiment Object
getUMICounts <- function(paths = NULL, sample.ids = NULL){
  
  message("1. Load Requirements")
  
  if(!require(SingleCellExperiment)){BiocManager::install("SingleCellExperiment")}
  suppressPackageStartupMessages(library(tidyverse))
  
  if(all(file.exists(paths))){
    message(paste0("All Provided Path Exists, proceed"))
  } else {
    stop(paste0("Following Path do not exists:", "\n", path[!file.exists(path)]))
  }
  
  # NB: Did not parallelize it. Few files. Should be fast enough as it is.
  message("2. Read Files")
  mycounts <- list()
  for(file in paths){
    mycounts[[file]] <- suppressMessages(read_delim(file, col_names = T, delim = "\t"))
  }
  mycounts <- purrr::reduce(mycounts, dplyr::full_join, by = c("gene")) %>%
    mutate_all(~replace_na(., 0)) %>%
    data.frame()
  
  row.names(mycounts) <- mycounts$gene
  mycounts$gene <- NULL
  colnames(mycounts) <- sample.ids
  
  message("3. Convert to SingleCellExperiment Object")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mycounts))
  
  return(sce)
  
}

################
# Median Genes # 
################

#' sceToMedianGeneCounts
#' 
#' Get the number of expressed gene per sample given one or more thresholds. Also returns the TPMs if rowRanges is present.
#' @param sce, singleCelObject
#' @param thresholds, vector of threshold to consider a gene as expressed
#' @param TPM, logical, should the TPM be calculated (will work only if RowRanges are present)
#' @param salmon, logical
sceToExpressedGenes <- function(sce = NULL, thresholds = c(0,1,5), salmon = FALSE, TPM = TRUE){
  
  suppressPackageStartupMessages(library(tidyverse))
  
  # 1. Get the median Gene count
  mycounts <- SummarizedExperiment::assay(sce)
  
  # 2. Count
  message("Get Median Counts")
  med.count <- lapply(thresholds, function(i) 
    data.frame(colSums(mycounts > i)) %>% 
      mutate(ID = colnames(mycounts))
  )
  
  # 3. TPM 
  # Check if rownames are similar to avoid dividing by the wrong width
  if(TPM == TRUE){
    message("Get Median TPM")
    if(salmon == FALSE){
      # Some of the matrix are too big. New approach to split them to avoid crash:
      gene.widths <- GenomicRanges::width(sce@rowRanges)/1000
      index <- split(1:ncol(mycounts), ceiling(seq_along(1:ncol(mycounts))/5000))
      
      med.gather <- list()
      for(i in 1:length(index)){
        print(paste0("Chunk: #", i))
        mycounts.rpk <- mycounts[,index[[i]]]
        totalReads <- colSums(mycounts.rpk)
        rpk <- mycounts.rpk/gene.widths
        # One value from 10_std.FS_C1_50000_2 is NA in one of the dowsampling (?). 
        rpk[is.na(rpk)] <- 0
        scaling <- colSums(rpk)/1e6
        tpm <- mapply("/", rpk, scaling)
        
        med.gather[[i]] <- lapply(thresholds, function(j) 
          data.frame(colSums(tpm > j, na.rm = T)) %>% 
            mutate(ID = colnames(mycounts.rpk))
        )
      }
      # Gather
      med.tpm <- list()
      for(j in 1:length(thresholds)){
        med.tpm[[j]] <- bind_rows(lapply(med.gather, function(x) x[[j]]))
      }
      
    } else if(salmon == TRUE){
      tpm <- SummarizedExperiment::assay(sce, "abundance")
      med.tpm <- lapply(thresholds, function(i) 
        data.frame(colSums(tpm > i)) %>% 
          mutate(ID = colnames(tpm)))
    }
  }
  # 4. Regroup data
  message("Regroup All")
  med.count <- purrr::reduce(med.count, left_join, by = "ID")
  
  if(TPM == TRUE){
    med.tpm <- purrr::reduce(med.tpm, left_join, by = "ID")
    med.all <- left_join(med.count, med.tpm, by = "ID") %>% dplyr::select(starts_with("ID"), everything())
    colnames(med.all) <- c("ID", paste0("count_", thresholds), paste0("TPM_", thresholds))
  } else {
    med.all <- med.count %>% dplyr::select(starts_with("ID"), everything())
    colnames(med.all) <- c("ID", paste0("count_", thresholds))
  }
  
  
  return(med.all)
}


#!/usr/bin/env Rscript

#' add_CB_CR_Barcodes
#' 
#' Add sequencing BC from index1/index2 (seq and qualities) to an unmapped bam
#' @param mybam, path to unmapped bam (FASTQtoSAM, Picard)
#' @param index1, path to index1.fastq.gz
#' @param index2, path to index2.fastq.gz

mybam="/home/vincent.hahaut/data_storage/transfer/single/121_1in5dT_96w_triton_C1.bam"
index1="/home/vincent.hahaut/data_storage/200904_NB551561_0030_AH5CW5BGXG/out/121_1in5-oligodT_FS_C1.I1.tmp.gz"
index2="/home/vincent.hahaut/data_storage/200904_NB551561_0030_AH5CW5BGXG/out/121_1in5-oligodT_FS_C1.I2.tmp.gz"


add_CB_CR_Barcodes <- function(mybam = NULL, index1 = NULL, index2 = NULL){
  
  require(ShortRead)
  require(Biostrings)
  require(tidyverse)
  
  # 1. Read Index FASTQ
  I1 <- ShortRead::readFastq(index1)
  I2 <- ShortRead::readFastq(index2)

  # 2. Index to Tibble
  indexes  <- data_frame(
    id = as.vector(ShortRead::id(I1)),
    seq = paste0(ShortRead::sread(I1), ShortRead::sread(I2)),
    qual = paste0(as.vector(I1@quality@quality), as.vector(I2@quality@quality))
  )
    
  # 3. Read BAM
  param <- ScanBamParam(what=c("qname","flag","rname", "strand", "pos", "qwidth", "mapq", "cigar","mrnm","mpos","isize","seq","qual"),tag=c("RG"))
  gr.aln <- readGAlignments(mybam, param=param, use.names=TRUE)
  
  
  ## setting up all bam flags (in my case this a BWA alignment and these are the default cols/flags)
  
  ## reading original bam file
  
  ## adding my custom bam tag "RX" (in this case you could use CB/UB etc..)
  ## in my case the UMI is part of the read name presented as "UM1_UM2_READNAME" you need to modify the regex below for your use case
  mcols(gr.aln)$RX <- sub("_.*","",sub("_","",mcols(gr.aln)$qname))
  
  ## this would be yours
  mcols(gr.aln)$CB <- sub("_.*","",mcols(gr.aln)$qname)
  mcols(gr.aln)$UB <- sub(".*_","",sub("#.*","",mcols(gr.aln)$qname))
  
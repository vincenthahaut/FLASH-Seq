#' FS_FASTQ_to_zUMIs_FASTQ
#' 
#' Smart-seq3-like processing is hardcoded in zUMIs pipeline. It requires R1 reads to harbor the SS3 adapter sequence at a defined position. 
#' Moreover, flexible search of UMI reads, spacer sequence and UMI in R2 are not available. Yet, their python script to reconstruct isoform sequences could be very useful to use.
#' This function takes a fastq file of FLASH-Seq data and transform it to a zUMIs compatible FASTQ. 
#' It includes a flexible search of the UMI in both R1/R2 reads and spacer trimming. Reads 1 of R2-UMI are used as reads 2 to restrand these reads.
#' In order to have zUMIs compatible data, the SS3 adapter is appended to every UMI reads, upstream of the UMI. 
#' 
#' @param spacer, 
#' @param ID, 
#' @param index1.path, 
#' @param index2.path, 
#' @param R1.path, 
#' @param R2.path, 

FS_FASTQ_to_zUMIs_FASTQ <- function(
  output.path = "/home/vincent.hahaut/Desktop/",
  spacer = "CTAAC",
  ID="316_DEEP",
  index1.path="/home/vincent.hahaut/data_storage/ORGANOIDS/out/316_FS_organoids_A2_S401_I1_001.fastq.gz",
  index2.path="/home/vincent.hahaut/data_storage/ORGANOIDS/out/316_FS_organoids_A2_S401_I2_001.fastq.gz",
  R1.path="/home/vincent.hahaut/data_storage/ORGANOIDS/out/316_FS_organoids_A2_S401_R1_001.fastq.gz",
  R2.path="/home/vincent.hahaut/data_storage/ORGANOIDS/out/316_FS_organoids_A2_S401_R2_001.fastq.gz"){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ShortRead))
  
  # 1. Read FASTQ Files
  message("1 - Read FASTQ")
  I1 <- readDNAStringSet(index1.path, format = "fastq", with.qualities = TRUE)
  I2 <- readDNAStringSet(index2.path, format = "fastq", with.qualities = TRUE)
  R1 <- readDNAStringSet(R1.path, format = "fastq", with.qualities = TRUE)
  R2 <- readDNAStringSet(R2.path, format = "fastq", with.qualities = TRUE)
  
  # 2. Tidy R1/R2 files
  message("2 - Tidy Information")
  reads <- data_frame(
    QNAME = str_split(names(R1), " ", simplify = TRUE)[,1],
    SEQ_R1 = as.vector(R1),
    QUAL_R1 = as.vector(mcols(R1)$qualities),
    SEQ_R2 = as.vector(R2),
    QUAL_R2 = as.vector(mcols(R2)$qualities)
  )
  
  # 3. UMI and spacer sequences
  
  # 3.1. Create the TSO pattern
  # = read starts with: adapter(trimmed) + UMI + SPACER + GG
  message("3 - Create Search Pattern")
  adapter_tso <- "AAGCAGTGGTATCAACGCAGAGT"
  pattern_adapter_R1 <- c(sapply(1:(nchar(adapter_tso)-2), function(x) str_sub(adapter_tso, nchar(adapter_tso)-x, nchar(adapter_tso)+1)), adapter_tso)
  pattern_adapter_R1 <- paste0("^", pattern_adapter_R1, "........", spacer, "GG")
  pattern_adapter_R1 <- paste0(pattern_adapter_R1, collapse = "|")
  
  pattern_adapter_R2 <- c("GAGT", "GT")
  pattern_adapter_R2 <- paste0("^", pattern_adapter_R2, "........", spacer, "GG")
  pattern_adapter_R2 <- paste0(pattern_adapter_R2, collapse = "|")
  
  identify_UMI_reads <- function(QNAME = NULL, SEQ = NULL, QUAL = NULL, pattern = NULL, TYPE = "R1"){
    
    # 3.2. Identify reads with UMI
    # No mismatch allowed
    # Trim sequences to remove FS adapter, UMI, spacer
    l_spacer_GG <- nchar(spacer) + 2
    
    # 3.2.1.1. Locate TSO position
    umi_position <- str_locate(pattern = pattern, string = SEQ)
    index_umi <- which(!is.na(umi_position[,1]))
    
    if(length(index_umi) > 1){
      # 3.2.1.2. Select UMI reads
      umi_position <- umi_position[index_umi,]
      read_width <- nchar(SEQ[index_umi])
      # 3.2.1.3. Trimmed read sequence
      seq_umi_trimmed <- str_sub(string = SEQ[index_umi], start = umi_position[,2],  end = read_width)
      qual_umi_trimmed <- str_sub(string = QUAL[index_umi], start = umi_position[,2],  end = read_width)
      # 3.2.1.4. UMI sequence
      umi_seq <- str_sub(string = SEQ[index_umi], start = umi_position[,2] - l_spacer_GG - 7, end = umi_position[,2] - l_spacer_GG)
      umi_qual <- str_sub(string = QUAL[index_umi], start = umi_position[,2] - l_spacer_GG - 7, end = umi_position[,2] - l_spacer_GG)
      
      data_frame(QNAME = QNAME[index_umi],
                 SEQ = seq_umi_trimmed,
                 QUAL = qual_umi_trimmed,
                 UB = umi_seq,
                 QU = umi_qual,
                 TYPE = TYPE) %>%
        return()
      
    } else {
      data_frame(QNAME = character(),
                 SEQ = character(),
                 QUAL = character(),
                 UB = character(),
                 QU = character(),
                 TYPE = TYPE) %>%
        return()
    }
    
  }
  
  R1_UMI_processed <- identify_UMI_reads(QNAME = reads$QNAME, SEQ = reads$SEQ_R1, QUAL = reads$QUAL_R1, pattern = pattern_adapter_R1, TYPE = "R1")
  R2_UMI_processed <- identify_UMI_reads(QNAME = reads$QNAME, SEQ = reads$SEQ_R2, QUAL = reads$QUAL_R2, pattern = pattern_adapter_R2, TYPE = "R2")
  
  message(paste0("R1-UMI: ", nrow(R1_UMI_processed)), " reads (", round(digits = 2, 100*nrow(R1_UMI_processed)/nrow(reads)), "%)")
  message(paste0("R2-UMI: ", nrow(R2_UMI_processed)), " reads (", round(digits = 2, 100*nrow(R2_UMI_processed)/nrow(reads)), "%)")
  
  if(nrow(R1_UMI_processed) > 0 | nrow(R2_UMI_processed) > 0){
    
    # 3.2.3. Combine R1 / R2-UMI
    message("6 - Combine UMI information")
    umi_reads <- bind_rows(R1_UMI_processed, R2_UMI_processed)
    
    # 4. Catch sequences with multiple UMI pattern detected
    # Quite rare (<0.001%) but duplicates read ids
    # Typically UMI pattern in both R1/R2 resulting from either R1 / R2 invasion (=short insert size) or retrotranscription errors
    message("7 - Clean-up UMI reads")
    duplicated_ids <- unique(umi_reads$QNAME[duplicated(umi_reads$QNAME)])
    umi_reads <- umi_reads[!umi_reads$QNAME %in% duplicated_ids,]
    
    # 5. Internal Reads
    # = Not duplicated IDs
    # = Not UMI reads
    message("8 - Get Internal Reads")
    internal_reads <- reads[! reads$QNAME %in% umi_reads$QNAME,]
    internal_reads <- internal_reads[! internal_reads$QNAME %in% duplicated_ids,]
    internal_reads$TYPE <- "INTERNAL"
    
    internal_reads <- dplyr::rename(internal_reads, SEQ = SEQ_R1,
                                    QUAL = QUAL_R1,
                                    MATE_SEQ = SEQ_R2,
                                    MATE_QUAL = QUAL_R2)
    # 6. Tidy UMI reads
    
    # 6.1. Reformat UMI read sequence to match zUMIs expectations 
    # zUMIs expects SS3 TSO format to classify a read as "UMI" and extract the sequence
    # Append the SS3 adapter / fake adapter quality
    message("9 - Add SS3 adapters to UMI reads")
    ss3_adapter_seq <- "ATTGCGCAATG"
    ss3_adapter_qual <- "EEEEEEEEEEE"
    
    UMI_cleaned <- mutate(umi_reads, 
                          SEQ = paste0(ss3_adapter_seq, UB, SEQ),
                          QUAL = paste0(ss3_adapter_seq, QU, QUAL))
    # 6.2. Add Mate
    # ==> Stranded R2-UMI
    # R2-UMI are sequenced from the opposite strand 
    # Exchange r1 and r2 reads from R2-UMI restrand R2-UMI reads
    # Doing this way keeps R2-UMI in "SEQ" which will be the R1_read later on
    UMI_cleaned <- mutate(UMI_cleaned, 
                          MATE_SEQ = ifelse(TYPE == "R1", reads$SEQ_R2[match(QNAME, reads$QNAME)], reads$SEQ_R1[match(QNAME, reads$QNAME)]),
                          MATE_QUAL = ifelse(TYPE == "R1", reads$QUAL_R2[match(QNAME, reads$QNAME)], reads$QUAL_R1[match(QNAME, reads$QNAME)]))
    
  } else {
    
    message("No UMI reads detected !!")
    
    UMI_cleaned <- data_frame()
    
    message("8 - Get Internal Reads")
    internal_reads <- reads
    internal_reads$TYPE <- "INTERNAL"
    
    internal_reads <- dplyr::rename(internal_reads, SEQ = SEQ_R1,
                                    QUAL = QUAL_R1,
                                    MATE_SEQ = SEQ_R2,
                                    MATE_QUAL = QUAL_R2)
    
  }
  
  # 7. Combine
  processed_reads <- bind_rows(internal_reads, UMI_cleaned)
  
  # 8. Add index information
  barcodes <- data_frame(QNAME = str_split(names(I1), " ", simplify = TRUE)[,1],
                         SEQ_I1 = as.vector(I1),
                         SEQ_I2 = as.vector(I2),
                         QUAL_I1 = as.vector(mcols(I1)$qualities),
                         QUAL_I2 = as.vector(mcols(I1)$qualities))
  
  processed_reads <- left_join(processed_reads, barcodes, by = "QNAME")
  
  # 9. Check-ups
  sum(duplicated(processed_reads$QNAME))
  nrow(reads) == nrow(processed_reads) + length(duplicated_ids)
  
  # 10. Write FASTQ file
  message("11 - Write Results")
  write.fastq <- function(seq = NULL, ids = NULL, qual = NULL, full.output.path = NULL){
    library(Biostrings)
    
    names(seq) <- ids
    dna.string <- DNAStringSet(seq)
    mcols(dna.string)$qualities <- BStringSet(qual)
    writeXStringSet(dna.string, full.output.path, compress = TRUE, format = "fastq")
    
  }
  
  write.fastq(seq = processed_reads$SEQ, 
              ids = processed_reads$QNAME, 
              qual = processed_reads$QUAL, 
              full.output.path = paste0(output.path, "/", ID, "_R1_zUMIs.fastq.gz"))
  
  write.fastq(seq = processed_reads$MATE_SEQ, 
              ids = processed_reads$QNAME, 
              qual = processed_reads$MATE_QUAL, 
              full.output.path = paste0(output.path, "/", ID, "_R2_zUMIs.fastq.gz"))
  
  write.fastq(seq = processed_reads$SEQ_I1, 
              ids = processed_reads$QNAME, 
              qual = processed_reads$QUAL_I1, 
              full.output.path = paste0(output.path, "/", ID, "_I1_zUMIs.fastq.gz"))
  
  write.fastq(seq = processed_reads$SEQ_I2, 
              ids = processed_reads$QNAME, 
              qual = processed_reads$QUAL_I2, 
              full.output.path = paste0(output.path, "/", ID, "_I2_zUMIs.fastq.gz"))
  
  message("12 - Finished!")
}


#### EXAMPLE: 

dir.create("/home/vincent.hahaut/data_storage/ORGANOIDS/FASTQ_zUMIs")

ID <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/out/", pattern = "R1_001.fastq.gz")
ID <- str_replace(string = ID, pattern = "_S.{1,3}_R1_001.fastq.gz$", "") 
R1.path <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/out/", pattern = "R1_001.fastq.gz", full.names = TRUE)
R2.path <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/out/", pattern = "R2_001.fastq.gz", full.names = TRUE)
I1.path <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/out/", pattern = "I1_001.fastq.gz", full.names = TRUE)
I2.path <- list.files("/home/vincent.hahaut/data_storage/ORGANOIDS/out/", pattern = "I2_001.fastq.gz", full.names = TRUE)

for(i in 1:length(ID)){
  
  message(paste0("ID: ", ID[i]))
  FS_FASTQ_to_zUMIs_FASTQ(
    output.path = "/home/vincent.hahaut/data_storage/ORGANOIDS/FASTQ_zUMIs/",
    spacer = "CTAAC",
    ID= ID[i],
    index1.path= I1.path[i],
    index2.path= I2.path[i],
    R1.path= R1.path[i],
    R2.path= R2.path[i])
  
}
  

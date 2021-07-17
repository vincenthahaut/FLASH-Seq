#' getDownstreamSeq
#' 
#' Deduplicate UMI reads, extract the sequence adjacent to the read start and add some custom annotations to summarize the invasion events.
#' @param mybam, /path/to/UMI_annotated.sorted.filtered.bam files
#' @param myid, sample ID
#' @param METHOD, FLASH-Seq or zUMIs depending on the pre-processing steps
#' @param genome, Loaded BSgenome file
getDownstreamSeq <- function(mybam = NULL, myid = NULL, METHOD = NULL, genome = genome){
  
  message(paste0("ID: ", myid))
  message("Read BAM File")
  
  # 1. Load Bam File
  # Allowed method:
  # zUMIs ==> Samples processed with zUMIs.
  # FLASH-Seq ==> Samples processed with umi_tools, STAR and featureCounts.
  if(METHOD == "zUMIs"){

    params <- ScanBamParam(tag = c("GE", "GI", "UB"), what=c("qname","rname","strand","pos","mapq", "cigar")) 
    bam <- scanBam(mybam, param = params)
    
    qwidth <- lapply(str_extract_all(pattern = "\\d*[MN]", bam[[1]]$cigar), function(x)
      str_replace_all(x, "[MN]", "") %>% 
        unlist() %>%
        as.numeric() %>%
        sum())
    
    gr <- makeGRangesFromDataFrame(
      data.frame(
        seqnames = bam[[1]]$rname,
        start = bam[[1]]$pos,
        end = bam[[1]]$pos + unlist(qwidth)-1,
        mapq = bam[[1]]$mapq,
        Assigned = "",
        exon = bam[[1]]$tag$GE,
        intron = bam[[1]]$tag$GI,
        cigar = bam[[1]]$cigar,
        UMI = bam[[1]]$tag$UB,
        strand = bam[[1]]$strand), 
      keep.extra.columns = TRUE) 
    
    # Keep only reads with a 8-bp UMI
    gr <- gr[nchar(as.character(gr$UMI)) == 8]
    
    gr$Assigned <- ifelse(!is.na(gr$exon) | !is.na(gr$intron), "Assigned", "Unassigned")
    
  } else if(METHOD == "FLASH-Seq"){
    
    params <- ScanBamParam(tag = c("XS", "XT"), what=c("qname","rname","strand","pos","qwidth","mapq", "cigar")) 
    bam <- scanBam(mybam, param = params)
    
    # Get the read end-of-mapping using CIGAR
    # While for read in +, start is enough to detect invasion, read in - requires the mapping end which is different from the BAM qwidth
    # Accounts for match (M) & skip (N = splicing)
    # Deletion (D), insertion (I) do not influence start/end
    # Soft-clipping are not mapping (S)
    qwidth <- lapply(str_extract_all(pattern = "\\d*[MN]", bam[[1]]$cigar), function(x)
      str_replace_all(x, "[MN]", "") %>% 
        unlist() %>%
        as.numeric() %>%
        sum())
    
    gr <- makeGRangesFromDataFrame(
      data.frame(
        seqnames = bam[[1]]$rname,
        start = bam[[1]]$pos,
        end = bam[[1]]$pos + unlist(qwidth)-1,
        mapq = bam[[1]]$mapq,
        UMI = str_split(bam[[1]]$qname, "_", simplify = TRUE)[,2],
        Assigned = bam[[1]]$tag$XS,
        Gene = bam[[1]]$tag$XT,
        strand = bam[[1]]$strand,
        cigar = bam[[1]]$cigar), 
      keep.extra.columns = TRUE)
  }
  
  # 2. Filter the GRanges
  # MAPQ > 5 / Main chromosomes
  gr <- gr[gr$mapq > 5]
  gr <- gr[seqnames(gr) %in% paste0("chr", 1:22)]
  
  if(length(gr) > 0){
    # 3. Get the sequence adjacent of the read start
    # 20bp before to account for potential mapping issues, long GGG's stretches created by the RT or longer spacer sequences.
    message("Get Upstream Sequences")
    
    gr$index <- 1:length(gr)
    
    gr.forward <- gr[strand(gr) == "+"]
    gr.forward.upstream <- GRanges(seqnames = seqnames(gr.forward), ranges = IRanges(start = start(gr.forward)-20, end = start(gr.forward)-1), strand = strand(gr.forward))
    gr.forward.upstream$index <- gr.forward$index

    gr.reverse <- gr[strand(gr) == "-"]
    gr.reverse.upstream <- GRanges(seqnames = seqnames(gr.reverse), ranges = IRanges(start = end(gr.reverse)+1, end = end(gr.reverse)+20), strand(gr.reverse))
    gr.reverse.upstream$index <- gr.reverse$index
    
    gr.upstream <- c(gr.forward.upstream, gr.reverse.upstream)
    gr.upstream <- gr.upstream[order(gr.upstream$index)]     # Reorder
    
    myseq.upstream <- Biostrings::getSeq(genome, gr.upstream)
    myseq.upstream <- as.data.frame(myseq.upstream)[,1]
    
    # 4. Combine results
    gr$upstream <- myseq.upstream
    gr.df <- mutate(as.data.frame(gr), 
                    READ_POS = paste0(seqnames, ":", start, "-", end),
                    READ_GROUP = paste0(seqnames, ":", start, "-", end, "-UMI:", UMI, "-upstream:", upstream)) 
    
    # 5. Collapse the data to remove duplicated
    # Based on the UMI, read mapping and adjacent sequence
    message("Collapse UMI positions")
    collapsed <- gr.df %>%
      group_by(READ_POS, upstream, UMI) %>%
      summarise(
        READ_STRAND = unique(strand),
        PCR_duplicates = n(),
        READ_GROUP = unique(READ_GROUP)) %>%
      ungroup() %>%
      mutate(ID = myid,
             UNIQUE_MOLECULES = n(),
             nREADS = length(gr))
    
    collapsed.dim <- nrow(collapsed)
    
    # 6. Add back the mapping info
    if(METHOD != "zUMIs"){
      collapsed <- left_join(collapsed, 
                             gr.df %>% 
                               select(READ_GROUP, Assigned, Gene) %>% distinct(), by = "READ_GROUP")
    } else {
      collapsed <- left_join(collapsed, 
                             gr.df %>% 
                               select(READ_GROUP, Assigned, exon, intron) %>% distinct(), by = "READ_GROUP")
    }
    
    # 7. In rare cases (<0.001%), reads with the same mapping positions (chr-start-end) can be assigned to different features
    # I think it comes from some secondary mappings not filtered (?). For now just exclude them.
    dubs <- unique(collapsed$READ_GROUP[duplicated(collapsed$READ_GROUP)])
    message(paste0("ID: ", myid, " - % Duplicated Multi-assigned Features: ", round(100*length(dubs) / collapsed.dim,6)))
    message(paste0("Filter ", length(collapsed$ID[collapsed$READ_GROUP %in% dubs]), " reads"))
    collapsed <- filter(collapsed, !READ_GROUP %in% dubs)
    
    # 8. Find if the sequence adjacent to the read matches the UMI
    message("Find Match")
    
    # 8.1. Create UMI pattern to search
    spacer <- case_when(grepl("5_SS3|6_SS3", myid) ~ "",
                        grepl("_AAGCA_", myid) ~ "AAGCA",
                        grepl("_CATCA_", myid) ~ "CATCA",
                        grepl("_CTGAC_", myid) ~ "CTGAC",                                    
                        grepl("_ATGAC_", myid) ~ "ATGAC",
                        grepl("_CTAAC_", myid) ~ "CTAAC",
                        grepl("SS3fwd", myid) ~ "",
                        grepl("STRToligoT_Spacer2", myid) ~ "ATAAC",
                        grepl("STRToligoT_Spacer1", myid) ~ "CAGCA",
                        grepl("STRToligoT_Fstso", myid) ~ "",
                        grepl("FSoligoT_Spacer2", myid) ~ "ATAAC",
                        grepl("FSoligoT_Spacer1", myid) ~ "CAGCA",
                        grepl("FSoligoT_FStso", myid) ~ "",
                        grepl("SS3oligoT_Spacer1", myid) ~ "CAGCA",
                        grepl("SS3oligoT_Fstso", myid) ~ "",
                        METHOD == "zUMIs" ~ "")
    
    collapsed$UMI.spacer <- paste0(collapsed$UMI, spacer, "GGG")
    
    # 8.2. Look for the pattern in the adjacent sequence
    # I observed several cases where there is not a perfect match. Likelihood of getting 5' mismatches is higher than internal mismatches in primers. 
    # But could be replaced later for partial matches (vmatchPattern)
    # Accounting for i mismatches (5')
    for(i in c(0,1,2,3)){
      
      umi <- as.character(substr(collapsed$UMI, 1+i, 8))
      umi.sp <- as.character(substr(collapsed$UMI.spacer, 1+i, unique(nchar(collapsed$UMI.spacer))))
      
      
      matches.umi <- mapply(function(x,y) grepl(x, y), umi, collapsed$upstream)
      matches.spa <- mapply(function(x,y) grepl(x, y), umi.sp, collapsed$upstream)
      #
      matches.umi.perc <- 100 * sum(matches.umi) / nrow(collapsed)
      matches.spa.perc <- 100 * sum(matches.spa) / nrow(collapsed)
      
      matches.df <- data.frame(matches.umi, matches.spa,matches.umi.perc, matches.spa.perc)
      mismatches <- ifelse(i != 0, paste0(".",i), "")
      
      colnames(matches.df) <- c(paste0("MATCH", mismatches), paste0("MATCH", mismatches, ".G"),
                                paste0("PERCENTAGE_MATCH", mismatches), paste0("PERCENTAGE_MATCH", mismatches, ".G"))
      
      print(paste0("Percentage Match (UMI) with ", i, " mismatches: ", round(100 * sum(matches.umi) / nrow(collapsed), 3)))
      print(paste0("Percentage Match (UMI-[Spacer]-GGG) with ", i, " mismatches: ", round(100 * sum(matches.spa) / nrow(collapsed), 3)))
      
      collapsed <- bind_cols(collapsed, matches.df)
      
    }
    
    # 9. Strand Invasion does not waited the UMI to occur. GGG is often enough
    motif <- substr(collapsed$upstream, 18,20)
    collapsed$GGG <- motif == "GGG"
    
    collapsed$PERCENTAGE_MATCH.GGG <- 100 * sum(motif == "GGG") / nrow(collapsed)
    
    # 10. Add Groups
    collapsed <-  mutate(collapsed, GROUP = case_when(grepl("5_SS3|6_SS3", ID) ~ "SS3",
                                                      grepl("_AAGCA_", ID) ~ "TSO-AAGCA\nSTRT-dT",
                                                      grepl("_CATCA_", ID) ~ "TSO-CATCA\nSTRT-dT",
                                                      grepl("_CTGAC_", ID) ~ "TSO-CTGAC\nSTRT-dT",
                                                      grepl("_ATGAC_", ID) ~ "TSO-ATGAC\nSTRT-dT",
                                                      grepl("_CTAAC_", ID) ~ "TSO-CTAAC\nSTRT-dT",
                                                      grepl("SS3fwd", ID) ~ "SS3\nHagemann-J.",
                                                      grepl("STRToligoT_Spacer2", ID) ~ "TSO-ATAAC\nSTRT-dT",
                                                      grepl("STRToligoT_Spacer1", ID) ~ "TSO-CAGCA\nSTRT-dT",
                                                      grepl("STRToligoT_Fstso", ID) ~ "TSO-FS-UMI\nSTRT-dT",
                                                      grepl("FSoligoT_Spacer2", ID) ~ "TSO-ATAAC\nFS-dT",
                                                      grepl("FSoligoT_Spacer1", ID) ~ "TSO-CAGCA\nFS-dT",
                                                      grepl("FSoligoT_FStso", ID) ~ "TSO-FS-UMI\nFS-dT",
                                                      grepl("SS3oligoT_Spacer1", ID) ~ "TSO-CAGCA\nSS3-dT",
                                                      grepl("SS3oligoT_Fstso", ID) ~ "TSO-FS-UMI\nSS3-dT",
                                                      METHOD == "zUMIs" ~ "zUMIs - SS3\nHagemann-J."))
  }
}


#### 0. RUN FUNCTION

# 1. Load the libraries
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(doParallel))

genome <- BSgenome.Hsapiens.UCSC.hg38

# 2. Prepare for multithreading
cl <- makeCluster(30)
registerDoSNOW(cl)
clusterCall(cl, function() library(tidyverse))
clusterCall(cl, function() library(Rsamtools))
clusterCall(cl, function() library(BSgenome.Hsapiens.UCSC.hg38))
  
# 3. Get paths
mypaths <- data_frame(
  bam.path = c("/home/vincent.hahaut/data_storage/SMART_SEQ3/HEK/zUMIs_output/demultiplexed/HEK.zUMIs.TTGTCGTGTACCACGA.demx.bam"),  # path/to/annotated/bam or zUMIs demultiplexed/*.demx.bam
  IDs =  c("TTGTCGTGTACCACGA"), # Associated sample IDs
  METHOD = c("zUMIs") # or FLASH-Seq
)

# 4. Process BAM files
res <- foreach(i = 1:nrow(mypaths)) %dopar% {
    getDownstreamSeq(mybam = mypaths$bam.path[i], myid = mypaths$IDs[i], METHOD = mypaths$METHOD[i], genome = genome)
  }
  
# 5. Stop the multi-threadhing
stopCluster(cl)
  
# 6. Combine results
res.combined <- bind_rows(res)
  
# 7. Save them
write_rds(res.combined, "~/bias.rds")


res.annotated.final <- read_rds("/home/vincent.hahaut/data_storage/FS_intermediate_file/bias.rds")

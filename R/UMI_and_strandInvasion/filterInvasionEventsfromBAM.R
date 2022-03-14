#!/usr/bin/env Rscript

############
# FUNCTION #
############

#' filterStrandInvasionToBAM
#' 
#' Get a list of read IDs containing UMI displaying a match to the read start adjacent sequence. Take this list and directly filter the bam file using RSamtools.
#' @param mybam, /path/to/UMI_annotated.sorted.filtered.bam files
#' @param myid, sample ID
#' @param METHOD, FLASH-Seq (=umi_tools) or zUMIs depending on the pre-processing steps (!!! zUMIs input were not thoroughly tested with this script)
#' @param fasta.genome, genome fasta file (need indexed fasta in the same folder ==> .fai)
#' @param minReads, minimum number of reads to consider the sample
#' @param mapq, minimum MAPQ to consider the read as UMI read
#' @param mismatches, maximum number of mismatches between the UMI and adjacent sequence
#' @param path.out, path/to/filteredID.txt (suffix only)

# mybam <- "/home/vincent.hahaut/data_storage/PBMC_PAIREDEND/STAR/382_FS-UMI_PBMC_A10/STAR/382_FS-UMI_PBMC_A10_ALL_Aligned.sortedByCoord.filtered.bam"
# myid <- "ID"
# genome <- "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"
# METHOD = "FLASH-Seq"
# minReads = 5000
# mismatches = 1
# mapq = 5
# path.out= "~/Desktop/"

filterStrandInvasionToBAM <- function(mybam = NULL, myid = NULL, METHOD = "FLASH-Seq", mismatches = 1, fasta.genome = genome, minReads = 5000, mapq = 5, path.out= "~/Desktop/"){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(Rsamtools))

  message("Read BAM File")
  
  # 1. Load Bam File
  # Allowed method:
  # zUMIs puts the UMI in the tags and umi_tools in the qname 
  # zUMIs ==> Samples processed with zUMIs (should work, not extensively tested on this particular function)
  # FLASH-Seq ==> Samples processed with umi_tools.
  message(paste0("Chosen Method: ", METHOD))
  if(METHOD == "zUMIs"){
    
    params <- ScanBamParam(tag = c("GE", "GI", "UB"), what=c("qname","rname","strand","pos","mapq", "cigar")) 
    bam <- scanBam(mybam, param = params)
    
    qwidth <- lapply(str_extract_all(pattern = "\\d*[MN]", bam[[1]]$cigar), function(x)
      str_replace_all(x, "[MN]", "") %>% 
        unlist() %>%
        as.numeric() %>%
        sum())
    
    # Saw one sample with no read at all associated to an exon tag (GE) - crashed the script
    if(any(sapply(bam[[1]]$tag, length) == 0)){
      bam[[1]]$tag <- lapply(bam[[1]]$tag, function(x) if(length(x) == 0) rep(NA, length(bam[[1]]$qname)) else x)
    }
    
    gr <- makeGRangesFromDataFrame(
      data.frame(
        qname = bam[[1]]$qname,
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
    
    gr$Assigned <- ifelse(!is.na(gr$exon) | !is.na(gr$intron), "Assigned", "Unassigned")
    
  } else if(METHOD == "FLASH-Seq"){
    
    params <- ScanBamParam(what=c("qname","rname","strand","pos","qwidth","mapq", "cigar")) 
    bam <- scanBam(mybam, param = params)
    
    # Get the read end-of-mapping using CIGAR
    # Accounts for match (M) & skip (N = splicing). Deletion (D), insertion (I) do not influence start/end. Soft-clipping are not mapping (S)
    qwidth <- lapply(str_extract_all(pattern = "\\d*[MN]", bam[[1]]$cigar), function(x)
      str_replace_all(x, "[MN]", "") %>% 
        unlist() %>%
        as.numeric() %>%
        sum())
    
    gr <- makeGRangesFromDataFrame(
      data.frame(
        qname = bam[[1]]$qname,
        seqnames = bam[[1]]$rname,
        start = bam[[1]]$pos,
        end = bam[[1]]$pos + unlist(qwidth)-1,
        mapq = bam[[1]]$mapq,
        UMI = str_split(bam[[1]]$qname, "_", simplify = TRUE)[,2],
        strand = bam[[1]]$strand,
        cigar = bam[[1]]$cigar), 
      keep.extra.columns = TRUE)
  
  }
  
  # 2. Filter the GRanges
  # MAPQ > 5 / Main chromosomes
  message(paste0("Total Reads: ", length(unique(gr$qname))))

  if(length(unique(gr$qname)) > minReads){
    
    # 3. Keep only UMI-reads
    gr <- gr[nchar(gr$UMI) > 0]
    message(paste0(myid, ": ", length(unique(gr$qname)), " UMI Reads"))
    
    # 4. Load Genome
    genome <- open(Rsamtools::FaFile(fasta.genome))
    
    # 5. Get the sequence adjacent of the read start
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
    
    # Addition of an adjacent sequence can create records bigger than the chromosome length = Crash
    # Only seen mitochondrial reads so far
    gr.index <- suppressMessages(read_tsv(col_names = FALSE, paste0(fasta.genome, ".fai")))
    gr.index <- makeGRangesFromDataFrame(gr.index, seqnames.field = "X1", start.field = "X2", end.field = "X2")
    start(gr.index) <- 1
    
    # Starts
    start(gr.upstream[start(gr.upstream) < 1]) <- 1 
    
    # Ends
    index <- findOverlaps(gr.upstream, gr.index, type = "within", ignore.strand = TRUE)
    correct.end <- gr.upstream[-queryHits(index)]
    chrom.end <- end(gr.index)[as.vector(match(seqnames(correct.end), seqnames(gr.index)))]
    end(gr.upstream[-queryHits(index)]) <- chrom.end
    
    gr.upstream <- gr.upstream[order(gr.upstream$index)]     # Reorder
    
    # Actual sequence extraction
    myseq.upstream <- Biostrings::getSeq(genome, gr.upstream)
    myseq.upstream <- as.data.frame(myseq.upstream)[,1]

    # 6. Combine results
    gr$upstream <- myseq.upstream

    gr.df <- mutate(as.data.frame(gr), 
                    READ_POS = paste0(seqnames, ":", start, "-", end),
                    READ_GROUP = paste0(seqnames, ":", start, "-", end, "-UMI:", UMI, "-upstream:", upstream)) 

    # 7. Compare UMI - adjacent sequence
    # max.distance controls the number of mismatches - not only 5' ones
    matches.umi <- mapply(function(x,y) agrepl(pattern = x, x = y, max.distance = as.numeric(mismatches)), as.character(gr.df$UMI), as.character(gr.df$upstream))

    message(paste0(myid, ": ", sum(matches.umi), " invasion events (", round(100*(sum(matches.umi)/length(unique(gr$qname))), 2), " % of the reads)"))
    
    # 8. Write a filtered IDs
    path.ids <- paste0(path.out, "/", myid, "_filteredInvasionID.txt")
    write.table(gr$qname[matches.umi], path.ids, sep= "\t", quote = F, row.names = F, col.names = FALSE)
  
    # 9. Write a filtered BAM
    # In case no index is found - might crash the filterBam()
    if(!file.exists(paste0(mybam, ".bai"))){indexBam(mybam)}
    
    message("Filter out invasion events")
    filterBam(mybam, 
              destination = paste0(path.out, "/", myid, "_filteredInvasion.bam"), 
              filter= FilterRules(list(subsetQNAME=function(x) {!x$qname %in% gr$qname[matches.umi]} )))
    write.table(data.frame(ID = myid, 
               Total_UMI_Reads = length(gr), 
               Invasion_events = sum(matches.umi),
               Invasion_events_perc = round(100*(sum(matches.umi)/length(gr)), 2)), 
               row.names = FALSE, sep = "\t", quote = FALSE,
               file = paste0(path.out, "/", myid, "_filteredInvasion.log.out"))
    message(paste0("Output Log: ", paste0(path.out, "/", myid, "_filteredInvasion.log.out")))
    message(paste0("Filter Bam: ", paste0(path.out, "/", myid, "_filteredInvasion.bam")))
  } else {
    message(paste0(myid, " not enough UMI reads - skip it"))
  }
  
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Stop - Requires 8 arguments - check the function for more info", call.=FALSE)
  message("arg1: path/to/bam")
  message("arg2: sample ID")
  message("arg3: method ('FLASH-Seq' or 'zUMIs')")
  message("arg4: number of mistmaches (recommended: 1)")
  message("arg5: path/to/genome.fa - genome used for the mapping")
  message("arg6: minimum number of total reads to process the sample (recommended: 10000)")
  message("arg7: minimum mapq to look for a match between UMI and adjacent sequence (recommended: 5)")
  message("arg8: path/to/output.bam")
} else if (length(args)==8) {
  # default output file
  mybam = args[1]
  myid = args[2]
  METHOD = args[3]
  mismatches = args[4]
  genome = args[5]
  minReads = args[6]
  mapq = args[7]
  path.out = args[8]
  
  message("####### ARGUMENTS ########")
  message(paste0("bam = ", args[1], " - Exists - ", file.exists(args[1])))
  message(paste0("Sample_ID = ", args[2]))
  message(paste0("Method = ", args[3]))
  message(paste0("Mistmatch = ", args[4]))
  message(paste0("genome = ", args[5]))
  message(paste0("min_reads = ", args[6]))
  message(paste0("mapq = ", args[7]))
  message(paste0("output_path = ", args[8]))
  message("##########################")
  
  
  filterStrandInvasionToBAM(mybam = mybam, 
                            myid = myid, 
                            METHOD = METHOD, 
                            mismatches = mismatches, 
                            genome = genome, 
                            minReads = minReads, 
                            mapq = mapq, 
                            path.out= path.out)
    
}

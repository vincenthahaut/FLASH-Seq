#!/usr/bin/R

# Extract Gene Length and GC content from a fasta file using a GTF as reference

if(!file.exists("/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.exon.GC_lengths.tsv")){
  # 99% from https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsamtools)
  
  GTFfile = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
  FASTAfile = "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"
  
  # 1. Load the annotation and reduce it
  GTF <- import.gff(GTFfile, format="gtf", genome="GRCh38", feature.type="exon")
  grl <- GenomicRanges::reduce(split(GTF, elementMetadata(GTF)$gene_name))
  reducedGTF <- unlist(grl, use.names=T)
  elementMetadata(reducedGTF)$gene_name <- rep(names(grl), elementNROWS(grl))
  
  # 2. Open the fasta file
  FASTA <- FaFile(FASTAfile)
  open(FASTA)
  
  # 3. Add the GC numbers
  elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
  elementMetadata(reducedGTF)$widths <- width(reducedGTF)
  
  # 4. Create a list of the ensembl_id/GC/length
  calc_GC_length <- function(x) {
      nGCs = sum(elementMetadata(x)$nGCs)
      width = sum(elementMetadata(x)$widths)
      c(width, nGCs/width)
  }
  output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_name), calc_GC_length))
  output <- rownames_to_column(as.data.frame(output), "gene_name")
  colnames(output) <- c("gene_name", "Length", "GC")
  
  write.table(output, file="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.exon.GC_lengths.tsv", quote = FALSE, sep="\t")
} else {
  
  gc.perc <- read.table(file="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.exon.GC_lengths.tsv", header = T, sep="\t")

}

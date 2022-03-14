# 0. Packages and FASTQ transcriptome
library(Biostrings)
library(ShortRead)
library(GenomicFeatures)
fa <- unlist(sread(readFasta("/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/TRANSCRIPTOME/gencode.v34.transcripts.fa.gz")))

## OPTION 1 - Transcriptome

# 1. Get the X-mer frequency
freq <- oligonucleotideFrequency(fa, width = 8, as.prob = TRUE)

View(as.data.frame(freq) %>% 
       rownames_to_column("oligo") %>% 
       mutate(SPACER = grepl(c("CAGCA|ATAAC|AAGCA|ATGAC|CATCA|CTAAC|CTGAC"), oligo)) %>%
       filter(!grepl("G$", oligo)) %>%
       #filter(grepl("GGG$", oligo), !grepl("GGGG$", oligo)) %>%
       arrange(freq) %>% 
       mutate(top = 1:n()))

## OPTION 1 - Genome
fa <- unlist(sread(readFasta("/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa")))

freqs <- list()
for(i in paste0("chr", 1:22)){
  freqs[[i]] <- oligonucleotideFrequency(genome[[i]], width = 8, as.prob = TRUE)
}
freq <- bind_rows(freqs)
freq <- apply(freq, 2 , function(x) mean(x))

View(as.data.frame(freq) %>% 
       rownames_to_column("oligo") %>% 
       mutate(SPACER = grepl(c("CAGCA|ATAAC|AAGCA|ATGAC|CATCA|CTAAC|CTGAC"), oligo)) %>%
       #filter(!grepl("G$", oligo)) %>%
       filter(grepl("GGG$", oligo), !grepl("GGGG$", oligo)) %>%
       arrange(freq) %>% 
       mutate(top = 1:n()))

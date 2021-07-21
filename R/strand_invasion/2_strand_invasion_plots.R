library(tidyverse)
library(ggseqlogo)

prefix <- "/path/to/output/"

# 0. Load the data generated using 1_strand_invasion_detection.R
res.annotated.final <- read_rds("bias.rds")

# 1. SeqLogo - 5-bp adjacent
res.annotated.final <- res.annotated.final %>%
  mutate(SUBGROUPS = case_when(GROUP == "SS3" ~ "SS3",
                               GROUP == "SS3\nHagemann-J." ~ "SS3\nHagemann-J.",
                               grepl("TSO-ATAAC|TSO-ATAAC|TSO-AAGCA|TSO-CATCA|TSO-CTGAC|TSO-ATGAC|TSO-CTAAC|TSO-CAGCA", GROUP) ~ "SPACER",
                               TRUE ~ "FS-UMI"))

p <- list()
for(i in unique(res.annotated.final$SUBGROUPS)){
  upstream <- filter(res.annotated.final, SUBGROUPS == i) %>%
    sample_n(size = 500000, replace = FALSE) 
  
  myseq <- substr(upstream$upstream, 15, 20)
  
  p[[i]] <- ggplot() + 
    ggseqlogo::geom_logo(seq_type = "DNA", 
                         method = "prob", 
                         col_scheme = "nucleotide", myseq) + 
    ggseqlogo::theme_logo(base_size = 24) 
  
}

ggsave(plot_grid(p[[1]], p[[2]],p[[3]], p[[4]], nrow = 1), filename = paste0(prefix, "seqlogo_upstream.GGG.tiff"), dpi= 300, height = 4, width = 20)

# 2. Intron / Exon 
# This part of the code is only compatible with data processed using the EXON/INTRON GTF generated with 0_Create_ExonIntron_Reference.R

# 2.1. Compute the % of deduplicated read per intron/exon accounting for the read/gene orientation
intron.exon.dat <- res.annotated.final %>% 
  filter(Assigned == "Assigned") %>%
  dplyr::select(GROUP, ID, READ_STRAND, Gene) %>%
  # The Gene column contain the following information:
  # EXON:Position-Gene_Orientation or INTRON:Position-Gene_Orientation
  mutate(TYPE = str_extract(Gene, "EXON|INTRON"),
         GENE_STRAND = str_extract(Gene, "\\+$|\\-$"),
         ANTISENSE = GENE_STRAND != READ_STRAND) %>%
  group_by(GROUP) %>%
  mutate(nREADS = n()) %>%
  group_by(GROUP, TYPE, ANTISENSE) %>%
  summarise(perc = 100*n()/unique(nREADS)) %>%
  ungroup() %>%
  # 
  mutate(GROUP = ifelse(!grepl("Hagemann", GROUP), str_replace(pattern = "\n", replacement = "+", GROUP), str_replace(pattern = "\n", replacement = " ", GROUP)),
         GROUP = factor(GROUP, levels = rev(c("TSO-AAGCA+STRT-dT","TSO-CATCA+STRT-dT","TSO-CTGAC+STRT-dT","TSO-ATGAC+STRT-dT","TSO-CTAAC+STRT-dT","TSO-ATAAC+STRT-dT","TSO-CAGCA+STRT-dT","TSO-ATAAC+FS-dT","TSO-CAGCA+FS-dT","TSO-CAGCA+SS3-dT","TSO-FS-UMI+SS3-dT","TSO-FS-UMI+STRT-dT","TSO-FS-UMI+FS-dT","zUMIs - SS3 Hagemann-J.","SS3","SS3 Hagemann-J."))),
         SUBGROUP = paste0(TYPE, " - ", ifelse(ANTISENSE == TRUE, "Discordant", "Concordant")),
         facet = case_when(GROUP %in% c("zUMIs - SS3 Hagemann-J.", "SS3", "SS3 Hagemann-J.") ~ "SS3",
                           grepl("TSO-ATAAC|TSO-ATAAC|TSO-AAGCA|TSO-CATCA|TSO-CTGAC|TSO-ATGAC|TSO-CTAAC|TSO-CAGCA", GROUP) ~ "SPACER",
                           TRUE ~ "FS-UMI")) %>%
  # From a check-up to verify if the invasion results are the same using zUMIs pipeline (==> They are)
  # Due to the reference zUMIs uses this data is not compatible with this plot
  filter(GROUP != "zUMIs - SS3 Hagemann-J.")


p.intron <- ggplot(intron.exon.dat, aes(x = GROUP, y = perc, fill = SUBGROUP)) + 
    geom_bar(stat = "identity", color = "black") +
    facet_grid(. ~facet, scales = "free", space = "free") + 
    theme_cowplot(font_size = 28) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          strip.background = element_blank(),
          axis.line = element_blank(),
          strip.text = element_text(face = "bold", size = 24)) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(breaks = seq(0,100,12.5)) + 
    xlab("Percentage 5'UMI-Reads") +
    ylab("") +
    guides(fill =guide_legend(nrow = 2))
  
ggsave(p.intron, filename = paste0(prefix, "intron_exon_UMI.tiff"), dpi= 300, height = 11, width = 13)

# 3. Strand Invasion %

# 3.1. Reformat the data for the graphics
dat <- res.annotated.final %>%
  filter(GROUP != "zUMIs - SS3\nHagemann-J.") %>% 
  select(PERCENTAGE_MATCH, PERCENTAGE_MATCH.1, PERCENTAGE_MATCH.2, PERCENTAGE_MATCH.3, PERCENTAGE_MATCH.G, PERCENTAGE_MATCH.1.G, PERCENTAGE_MATCH.2.G, PERCENTAGE_MATCH.3.G, ID, GROUP) %>%
  distinct() %>%
  gather(mismatch, percentage, -GROUP, -ID) %>%
  mutate(SUBGROUPS = case_when(GROUP %in% c("zUMIs - SS3\nHagemann-J.", "SS3", "SS3\nHagemann-J.") ~ "SS3",
                               grepl("TSO-ATAAC|TSO-ATAAC|TSO-AAGCA|TSO-CATCA|TSO-CTGAC|TSO-ATGAC|TSO-CTAAC|TSO-CAGCA", GROUP) ~ "SPACER",
                               TRUE ~ "FS-UMI"),
         mismatchN = case_when(mismatch %in% c("PERCENTAGE_MATCH", "PERCENTAGE_MATCH.G") ~ 0,
                               mismatch %in% c("PERCENTAGE_MATCH.1", "PERCENTAGE_MATCH.1.G")  ~ 1,
                               mismatch %in% c("PERCENTAGE_MATCH.2", "PERCENTAGE_MATCH.2.G")  ~ 2,
                               mismatch %in% c("PERCENTAGE_MATCH.3", "PERCENTAGE_MATCH.3.G")  ~ 3),
         GGG = case_when(mismatch %in% c("PERCENTAGE_MATCH.G", "PERCENTAGE_MATCH.1.G", "PERCENTAGE_MATCH.2.G", "PERCENTAGE_MATCH.3.G") & SUBGROUPS == "SPACER" ~ "5'-UMI-(Spacer)-GGG-3'",
                         mismatch %in% c("PERCENTAGE_MATCH.G", "PERCENTAGE_MATCH.1.G", "PERCENTAGE_MATCH.2.G", "PERCENTAGE_MATCH.3.G") & SUBGROUPS %in% c("SS3", "FS-UMI") ~ "5'-UMI-GGG-3'",
                         TRUE ~ "5'-UMI-3'"),
         facet = factor(paste0(GGG, "\n", SUBGROUPS), levels = c(
           "5'-UMI-3'\nSS3", "5'-UMI-3'\nFS-UMI","5'-UMI-3'\nSPACER",
           "5'-UMI-GGG-3'\nSS3",  "5'-UMI-GGG-3'\nFS-UMI", "5'-UMI-(Spacer)-GGG-3'\nSPACER" )),
         percentage = as.numeric(percentage),
         SUBGROUPS = factor(SUBGROUPS, levels = c("SS3", "FS-UMI", "SPACER"))) %>%
  # Mean / Sd per mismatch, GROUP
  group_by(mismatchN, GGG, GROUP, SUBGROUPS, facet) %>%
  summarise(m = mean(percentage),
            s = sd(percentage))

# 3.2.1. Including Spacer/GGG
p.all.GGG <- filter(dat, facet %in% c("5'-UMI-GGG-3'\nSS3",  "5'-UMI-GGG-3'\nFS-UMI", "5'-UMI-(Spacer)-GGG-3'\nSPACER")) %>%
  ggplot(aes(x = mismatchN, ymin = m-s, ymax = m+s, y = m, color = GROUP)) +
  geom_pointrange() +
  geom_line(aes(y = m, group = GROUP)) +
  facet_wrap(~SUBGROUPS, nrow = 1) +
  scale_color_manual(values = c("SS3" = "#e6194B",
                                "TSO-ATAAC\nSTRT-dT"= "#3cb44b",
                                "TSO-CAGCA\nSTRT-dT"= "#ffe119",
                                "TSO-FS-UMI\nSTRT-dT"= "#4363d8",
                                "TSO-ATAAC\nFS-dT"= "#f58231",
                                "TSO-CAGCA\nFS-dT"= "#42d4f4",
                                "TSO-FS-UMI\nFS-dT"= "#f032e6",
                                "TSO-CAGCA\nSS3-dT"= "#fabed4",
                                "TSO-FS-UMI\nSS3-dT"= "#469990",
                                "SS3\nHagemann-J." = "#dcbeff",
                                "zUMIs - SS3\nHagemann-J." = "#B3B3B3",
                                "TSO-AAGCA\nSTRT-dT" = '#9A6324',
                                "TSO-CATCA\nSTRT-dT" ='#bfef45',
                                "TSO-CTGAC\nSTRT-dT" = '#800000',
                                "TSO-ATGAC\nSTRT-dT" ='#000075',
                                "TSO-CTAAC\nSTRT-dT" = '#808000')) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(0,3)) +
  xlab("Number of Consecutive 5' Mismatches in UMI") + 
  ylab("% Upstream Sequence Match") 

p.all <- filter(dat, facet %in% c("5'-UMI-3'\nSS3", "5'-UMI-3'\nFS-UMI","5'-UMI-3'\nSPACER")) %>%
  ggplot(aes(x = mismatchN, ymin = m-s, ymax = m+s, y = m, color = GROUP)) +
  geom_pointrange() +
  #ylim(c(0,20.5)) +
  geom_line(aes(y = m, group = GROUP)) +
  facet_wrap(~SUBGROUPS, nrow = 1) +
  scale_color_manual(values = c("SS3" = "#e6194B",
                                "TSO-ATAAC\nSTRT-dT"= "#3cb44b",
                                "TSO-CAGCA\nSTRT-dT"= "#ffe119",
                                "TSO-FS-UMI\nSTRT-dT"= "#4363d8",
                                "TSO-ATAAC\nFS-dT"= "#f58231",
                                "TSO-CAGCA\nFS-dT"= "#42d4f4",
                                "TSO-FS-UMI\nFS-dT"= "#f032e6",
                                "TSO-CAGCA\nSS3-dT"= "#fabed4",
                                "TSO-FS-UMI\nSS3-dT"= "#469990",
                                "SS3\nHagemann-J." = "#dcbeff",
                                "zUMIs - SS3\nHagemann-J." = "#dcbeff",
                                "TSO-AAGCA\nSTRT-dT" = '#9A6324',
                                "TSO-CATCA\nSTRT-dT" ='#bfef45',
                                "TSO-CTGAC\nSTRT-dT" = '#800000',
                                "TSO-ATGAC\nSTRT-dT" ='#000075',
                                "TSO-CTAAC\nSTRT-dT" = '#808000',
                                "LA-TSO-CAGCA\nSTRT-dT" = "#dcbeff")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,1,2,3), limits = c(0,3)) +
  xlab("Number of Consecutive 5' Mismatches in UMI") + 
  ylab("% Upstream Sequence Match") 

ggsave(p.all, filename = paste0(prefix, "StrandInvasion_Mismatch.tiff"), dpi= 300, width = 14, height = 5.5)
ggsave(p.all.GGG, filename = paste0(prefix, "StrandInvasion_Mismatch.GGG.tiff"), dpi= 300, width = 14, height = 5.5)
  




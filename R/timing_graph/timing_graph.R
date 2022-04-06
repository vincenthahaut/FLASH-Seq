library(gggenes)  # Arrow charts are not well covered by ggplot but worked well with this one.
library(tidyverse)
library(cowplot)

timing <- read_tsv("/home/vincent.hahaut/Desktop/FLASH-Seq/R/timing_graph/FS_timing.txt")

# 0. Should you produce the first main main figure or the Low Amplification comparison ? 
LOWAMP = "ALL"

# 1. Tidy the data
mytimes <- timing %>%
  mutate(Step = factor(Step, levels = c("Start", "Denature", "RT", "RT-PCR", "PCR", "Clean-up", "QC", "Dilution", "Tagmentation", "Library Pooling", "Library Amplification", "Pooling and Clean-up"))) %>%
  filter(Step != "Start") %>%
  group_by(METHOD) %>%
  mutate(strand = "forward", 
         start = cumsum(duration)-duration,
         start = ifelse(start-10 > 0, start - 10, 0),
         end = cumsum(duration),
         totalTime = sum(duration)) %>%
  ungroup() %>%
  arrange(METHOD, desc(start)) %>%
  # The for loop below mess up the factor levels, easier to do it like this:
  mutate(METHOD = case_when(
    METHOD == "FS-LA" ~ "1-FS-LA",
    METHOD == "FS" ~ "2-FS",
    METHOD == "SS3" ~ "3-SS3",
    METHOD == "SS2" ~ "4-SS2",
    METHOD == "SSsc" ~ "5-SSsc"))

# 2. Add the x-axis
if(LOWAMP == FALSE){
  mytimes.tpm <- filter(mytimes, METHOD != "1-FS-LA")
  pA <- ggplot(mytimes.tpm) +
    geom_vline(xintercept = seq(0,600,60), linetype = "dashed", color = "darkgrey") +
    scale_x_continuous(breaks = seq(0,600,60), labels = 0:10, "Hours") 
} else if(LOWAMP == TRUE) {
  mytimes.tpm <- filter(mytimes, METHOD %in% c("1-FS-LA", "2-FS"))
  pA <- ggplot(mytimes.tpm) +
    geom_vline(xintercept = seq(0,420,60), linetype = "dashed", color = "darkgrey") +
    scale_x_continuous(breaks = seq(0,420,60), labels = 0:7, "Hours") 
} else if(LOWAMP == "ALL"){
  mytimes.tpm <- mytimes
  pA <- ggplot(mytimes.tpm) +
    geom_vline(xintercept = seq(0,600,60), linetype = "dashed", color = "darkgrey") +
    scale_x_continuous(breaks = seq(0,600,60), labels = 0:10, "Hours") 
  
}

# 3. Create the theme of the graph
pA <- pA +theme_cowplot() +
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 28),
        legend.position = "none",
        axis.text.x = element_text(size = 22), 
        axis.text.y = element_text(size = 22, face = "bold"), 
        legend.title = element_blank()) +
  ylab("") +
  scale_fill_manual(values = 
                      c("Denature" = '#1B9E77',  
                        "RT" = "#D95F02", 
                        "RT-PCR" = "#7570B3", 
                        "PCR" = "#E7298A", 
                        "Clean-up" = "#66A61E", 
                        "QC" = "#E6AB02", 
                        "Dilution" = "#A6761D", 
                        "Tagmentation" = "#666666", 
                        "Library Pooling" = "#E41A1C", 
                        "Library Amplification" = "#377EB8", 
                        "Pooling and Clean-up" = "#FF7F00")) 

# 4. Add the blocks one-by-one as arrows using gggenes
# The loop allows the addition of each step on top of the previous, creating the white separations
for(row in 1:nrow(mytimes.tpm)){
  pA <- pA +
    geom_gene_arrow(data = mytimes.tpm[row,], aes(xmin = start, xmax = end, y = METHOD, fill = Step), 
                    color = "white", 
                    arrow_body_height = unit(6, "mm"), 
                    arrowhead_height = unit(6, "mm"), 
                    arrowhead_width = unit(4, "mm"))
  
}

# 5. The first 3 min-Denaturation does not appear well (too small) - putting it on top helps
pA <- pA + geom_gene_arrow(data = mytimes.tpm[mytimes.tpm$METHOD == "Denature",], aes(xmin = start, xmax = end, y = METHOD, fill = Step),
                           color = "white", 
                           arrow_body_height = unit(6, "mm"), 
                           arrowhead_height = unit(6, "mm"), 
                           arrowhead_width = unit(4, "mm"))

# 6. Change the y-label
# 7. Save the results
if(LOWAMP == FALSE){
  ggsave(pA + 
           scale_y_discrete(labels=c("FS", "SS3", "SS2", "SSsc")) + 
           theme(legend.position = "bottom") + 
           guides(colour = guide_legend(nrow = 3, override.aes = list(size=12, stroke = 2), title = "Start")),
         filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/timing.tiff", dpi = 450, width = 14, height = 6, bg = "white")
} else if(LOWAMP == TRUE){
  ggsave(pA + 
           scale_y_discrete(labels=c("FS-LA", "FS")) + 
           guides(colour = guide_legend(override.aes = list(size=12, stroke = 2), title = "Start")), 
         filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/timing_FSonly.tiff", dpi = 450, width = 12, height = 3, bg = "white")
} else {
  ggsave(pA + 
           scale_y_discrete(labels=rev(c("SSsc", "SS2", "SS3", "FS", "FS-LA"))) + 
           guides(colour = guide_legend(override.aes = list(size=12, stroke = 2), title = "Start")), 
         filename = "/home/vincent.hahaut/data_storage/FS_intermediate_file/FIGURES/timing_FSonly.ALL.pdf", dpi = 450, width = 12, height = 4.5, bg = "white")
  
}

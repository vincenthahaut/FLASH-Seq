library(tidyverse)
library(kableExtra)
library(knitr)

conditions <- read_tsv("/home/vincent.hahaut/Desktop/FLASH-Seq/R/extended_files/1_summary_table_conditions.txt")
mycols <- conditions$Effect
mycols <- ifelse(mycols == "grey70", "grey", mycols)
mycols <- ifelse(mycols == "goldenrod1", "yellow", mycols)
mycols <- ifelse(mycols == "darkolivegreen1", "lightgreen", mycols)

conditions$Effect <- rep("", nrow(conditions))



condition1 <- filter(conditions, `Protocol step` %in% c("Lysis Mix", "Lysis & mRNA denaturation", "RT-PCR reaction"))
  condition1 %>%
  dplyr::select(-`Protocol step`) %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  row_spec(0, font_size = 20, bold = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), stripe_color = "black", font_size = 14) %>%
  column_spec(2, color = "black",
              background = mycols) %>%
  pack_rows("Lysis Mix", min(which(condition1$`Protocol step` == "Lysis Mix")), max(which(condition1$`Protocol step` == "Lysis Mix"))) %>%
  pack_rows("Lysis & mRNA denaturation", min(which(condition1$`Protocol step` == "Lysis & mRNA denaturation")), max(which(condition1$`Protocol step` == "Lysis & mRNA denaturation"))) %>%
  pack_rows("RT-PCR reaction", min(which(condition1$`Protocol step` == "RT-PCR reaction")), max(which(condition1$`Protocol step` == "RT-PCR reaction"))) %>%
  footnote(
    general = "If not clearly stated, concentration of reagents and reaction conditions are the same as in the standard FLASH-seq protocol. ") %>%
  save_kable("/home/vincent.hahaut/Desktop/table1.png", density = 600)

condition2 <- filter(conditions, `Protocol step` %in% c("RT-PCR Mix"))
  
condition2 %>%
  filter(`Protocol step` %in% c("RT-PCR Mix")) %>%
  dplyr::select(-`Protocol step`) %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Arial") %>%
  row_spec(0, font_size = 20, bold = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), stripe_color = "black", font_size = 14) %>%
  column_spec(2, color = "black",
              background = mycols) %>%
  pack_rows("RT-PCR Mix", min(which(condition2$`Protocol step` == "RT-PCR Mix")), max(which(condition2$`Protocol step` == "RT-PCR Mix"))) %>%
  footnote(number = c(
    "Concentration in the Lysis Mix; ", 
    "As FLASH-seq is performed in the KAPA HiFi mix and not in a standard RT mix, we could only add NaCl on top and not instead of the reaction buffer.
    ; ",
    "Even if the TSO is added with the Lysis Mix, the value reported here refers to the final concentration in the RT-PCR reaction. ; ", 
    "Due to volume constraints, Ficoll addition required the removal of betaine from the RT-PCR reaction mix. ; ", 
    "KAPA HiFi HotStart readyMix contains 0.3 mM dNTPs. FS Lysis Mix contains 1.2 mM dNTPs. ; "),
    general = "If not clearly stated, concentration of reagents and reaction conditions are the same as in the standard FLASH-seq protocol. ") %>%
  save_kable("/home/vincent.hahaut/Desktop/table2.png", density = 600)
# Strand Invasion Detection

## The procedure used to explore strand-invasion events contains three steps:

1. Create a GTF file containing the collapsed exon-intron positions of protein-coding genes using **0_Create_ExonIntron_Reference.R**.

2. Process the FASTQ files to extract UMIs, map and assign to a feature the 5' UMI reads. See **4_SINGLE_END_UMI.sh** for examples.

3. Parse the annotated BAM files uisng **1_strand_invasion_detection.R** to deduplicate UMI reads and extract the sequence adjacent to the read start.

## Notes regarding 1_strand_invasion_detection.R:

* This script is not designed to remove putative invasion events. However, it will provide an overview of such events.

* It is currently quite slow as it contains a lot of unecessary check-ups which were placed to explore the results. 

* Two inputs are accepted:
  - BAM files processed with **4_SINGLE_END_UMI.sh**.
  - zUMIs BAM files (https://github.com/sdparekh/zUMIs). 

Please note that outputs using zUMIs BAMs will not be compatible with my scripts to look for discordant/concordant reads as they have been designed around the GTF created in **0_Create_ExonIntron_Reference.R** but zUMIs creates a new GTF on-the-fly.

* It should work on paired-end reads but was not thoroughly tested. As the UMI are in theory located in R1, please use only single-read inputs for now. 

* BAM files have been pre-filtered using *samtools view -F 260*. 

* The script contains some hardcoded values for the human genome (e.g., chromosome number) and FLASH-Seq (e.g., spacer, ...). 




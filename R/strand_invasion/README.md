# Strand Invasion Detection

## Process

1. Create a GTF file containing the collapsed exon-intron positions of protein-coding genes using **0_Create_ExonIntron_Reference.R**.

2. Process the FASTQ files to extract UMIs, map and assign to a feature the 5' UMI reads. See **4_SINGLE_END_UMI.sh** for examples or use the polyvalent zUMIs (https://github.com/sdparekh/zUMIs).

3. Parse the annotated BAM files uisng **1_strand_invasion_detection.R** to deduplicate UMI reads and extract the sequence adjacent to the read start.


## Notes:

**GOAL:** This script is not designed to remove putative invasion events. However, it will provide an overview of the affected reads.

* Two inputs are accepted:
  - BAM files processed with **4_SINGLE_END_UMI.sh**.
  - zUMIs BAM files (https://github.com/sdparekh/zUMIs). 

Please note that zUMIs output BAMs are not compatible with the function used to look for protein-coding discordant/concordant reads as it was designed to work with samples processed using the GTF created in **0_Create_ExonIntron_Reference.R**. zUMIs creates a new GTF on-the-fly.

* It should work on paired-end reads but was not thoroughly tested. As UMIs are in theory located in R1, please use only single-read inputs for now. 

* BAM files have been pre-filtered using *samtools view -F 260*. 

* The script contains some hardcoded values for this study [human genome (e.g., chromosome number) and FLASH-Seq groups (e.g., spacer, ...)]

* Two bam files produced using zUMIs (HEK and HCA) are provided as example.

* I am aware that 1_strand_invasion_detection.R is currently quite slow. Improvments can still be made but in the meanwhile please parallelise the function as described in the script. 



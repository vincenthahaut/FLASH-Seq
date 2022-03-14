# Ongoing-development

## FS_FASTQ_to_zUMIs_FASTQ

**==> Not extensively tested !!!**

Smart-seq3 processing is hardcoded in zUMIs pipeline. It requires R1 reads to harbor the SS3 adapter sequence at a defined position. In addition, zUMIs only supports a fixed UMI position in R1 (and not R2) or FS spacer sequence. Yet, their python script to reconstruct isoform sequences could be very useful to use. This function takes demultiplexed FASTQ files (R1/R2/I1/I2) and creates a "fake" zUMIs compatible FASTQ:

1. Identify UMI reads (R1/R2)
2. Extract UMI and remove spacer sequence
3. Append SS3 adapter sequence upstream of the UMI
4. Reconciliate R1/R2 UMI reads

The zUMIs-comaptible FASTQ files can then be all concatenated for zUMIs. 

#!/bin/bash

# 0. Define paths
BASECALL_DIR="/home/vincent.hahaut/data_storage/220304_NB551561_0086_AHKCVLBGXL/Data/Intensities/BaseCalls/"
OUTPUT_DIR="/home/vincent.hahaut/data_storage/220304_NB551561_0086_AHKCVLBGXL/out/"
SAMPLESHEET="/home/vincent.hahaut/Desktop/GITHUB/IOB/sampleSheet/Pool94_04032022_FSimprovment_AHKCVLBGXL.csv"

# 1. Run BCL2FASTQ
# Keep index and regroup lane info
ulimit -n 10000
bcl2fastq --input-dir $BASECALL_DIR --output-dir $OUTPUT_DIR --sample-sheet $SAMPLESHEET --create-fastq-for-index-reads --no-lane-splitting

# 2. Sanity Check: What are the indexes left ? 
zcat out/Undetermined_S0_I1_001.fastq.gz | awk -F' 1:N:0:' 'NR%4==1{print $2}' | sort | uniq -c > left_index.txt
sort -k1,1 left_index.txt

# 3. How many reads per samples ? 
# Divide values by 4
for file in ./out/*R1*
do
zcat $file | wc -l 
done


#### Quick detection of UMI proportion
for file in FS*R1*; do echo $file; zcat $file | wc -l ; zgrep -c GT........CTAACGG $file; done
for file in SS3*R1*; do echo $file; zcat $file | wc -l ; zgrep -c TG........CGTACGG $file; done

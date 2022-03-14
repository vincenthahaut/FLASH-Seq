#!/bin/bash

# 0. Define paths
OUTPUTREF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/STAR/"
#OUTPUTREF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/STAR_100bp/"
FASTA="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
STAR="/home/vincent.hahaut/binaries/STAR-2.7.3a/bin/Linux_x86_64/STAR"

# 1. Generate the STAR reference genome
mkdir $OUTPUTREF
cd $OUTPUTREF

scp $FASTA .
scp $GTF .

$STAR --runThreadN 30 \
--runMode genomeGenerate \
--genomeDir $OUTPUTREF \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 74 # Read Length - 1
# --sjdbOverhang 29 # Read Length - 1
# --sjdbOverhang 99 # Read Length - 1
# --sjdbOverhang 49 # Read Length - 1


# Salmon indexes where obtained from:
# http://refgenomes.databio.org/v3/assets/splash/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_sa_index?tag=default
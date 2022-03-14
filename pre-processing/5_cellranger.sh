#!/bin/bash

# 1. Index Genome
mkdir /home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/hg38_gencode_10x/
cd /home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/hg38_gencode_10x/
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
FILTERED=hg38_gencode_FILTERED.gtf
FASTA="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"
GENOME=hg38_gencode

# 2. Filter the GTF to keep only the genes of interest
export PATH=/home/vincent.hahaut/data_storage/binaries/cellranger-6.1.1:$PATH

cellranger mkgtf $GTF $FILTERED \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

# 3. Create the reference
cellranger mkref \
--genome=$GENOME \
--fasta=$FASTA \
--genes=$FILTERED

mv $GENOME ../


# 4. Demultiplex the sequencing results
export PATH=/home/vincent.hahaut/data_storage/binaries/cellranger-6.1.1:$PATH
SAMPLESHEET="/home/vincent.hahaut/Desktop/GITHUB/SampleSheet/211105_NB551561_0068_AHCM3WAFX3_demultiplexing_sampleSheet.txt"
ILLUMINA_PATH="/home/vincent.hahaut/data_storage/ORGANOIDS_W18/FASTQ/211105_NB551561_0068_AHCM3WAFX3"
OUT="/home/vincent.hahaut/data_storage/ORGANOIDS_W18/FASTQ/211105_NB551561_0068_AHCM3WAFX3/out/"

mkdir $OUT
cd $OUT

# --qc Compute 10x and sequencing metrics
# --run BCL Illumina Folder
# --id Folder created by cellranger-atac
$cellranger mkfastq --run $ILLUMINA_PATH --id 211105_NB551561_0068_AHCM3WAFX3 --csv $SAMPLESHEET

# 5. Run CellRanger
# The retinal organoid sample has been been sequenced twice.
# Combine the results in cellranger count
export PATH=/home/vincent.hahaut/data_storage/binaries/cellranger-6.1.1:$PATH
cellranger count --localmem=50 --fastqs="/home/vincent.hahaut/data_storage/ORGANOIDS_W18/FASTQ/BSSE_QGF_187185_HL7Y2DRXY,/home/vincent.hahaut/data_storage/ORGANOIDS_W18/FASTQ/211105_NB551561_0068_AHCM3WAFX3" --sample="BSSE_QGF_187185_HL7Y2DRXY_2_Pool53_508_SI-GA-C10,BSSE_QGF_187185_HL7Y2DRXY_1_Pool53_508_SI-GA-C10" --id="ORGANOIDS_W18" --transcriptome="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/hg38_gencode/" --expect-cells=8000

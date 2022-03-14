#!/bin/bash

# 0. Paths
REFERENCE="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"
STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR_100bp_SJ/"
PICARD="/home/vincent.hahaut/data_storage/binaries/PICARD/picard.jar"   
GATK="/home/vincent.hahaut/data_storage/binaries/gatk-4.2.2.0/gatk"
cd "/home/vincent.hahaut/data_storage/ORGANOIDS"


# 1. STAR reference
# Create a two-pass reference genome containing selected splice-junctions from the first pass
# Run only one time

# SJ=`ls /home/vincent.hahaut/data_storage/ORGANOIDS/STAR/*/STAR/*ALL*_SJ.out.tab`
# GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
# 0. Clean up the SJ file
# See https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I
# cd $STAR_REF
# cat $SJ > SJ.all
# cat SJ.all | grep -v chrM | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ.filtered

# STAR --runMode genomeGenerate --genomeDir "$STAR_REF" --genomeFastaFiles $REFERENCE --sjdbGTFfile $GTF --runThreadN 30 --sjdbOverhang 99 --sjdbFileChrStartEnd SJ.filtered

STAR --genomeLoad LoadAndExit --genomeDir $STAR_REF

ulimit -n 100000

mkdir "/home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS"

# 2. Variant Calling
cat "/home/vincent.hahaut/data_storage/ORGANOIDS/xah" |
while IFS=$'\t' read -r -a myArray; do

ID="${myArray[0]}"
echo "$ID"
MYFASTQ_R1=/home/vincent.hahaut/data_storage/ORGANOIDS/STAR/"$ID"/FASTQ/allreads.R1.fq.gz
MYFASTQ_R2=/home/vincent.hahaut/data_storage/ORGANOIDS/STAR/"$ID"/FASTQ/allreads.R2.fq.gz

mkdir /home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/"$ID"
cd /home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/"$ID"

for DOWN in 5000 10000 20000 40000 50000 75000 100000 125000 250000 375000 500000 750000
do
      echo "$DOWN"
      if [ `zcat "$MYFASTQ_R1" | grep -c ^@`  -gt "$DOWN" ]
      then

            # 2.0. Resampling
            seqtk sample -s42 "$MYFASTQ_R1" "$DOWN" > R1."$DOWN".fastq
            seqtk sample -s42 "$MYFASTQ_R2" "$DOWN" > R2."$DOWN".fastq

            java -jar $PICARD FastqToSam \
                  F1=R1."$DOWN".fastq \
                  F2=R2."$DOWN".fastq \
                  O=unaligned_reads.bam \
                  SM="$ID" \
                  PL="Illumina" \
                  RG="$ID"

            # 2.1. Re-map 2-pass using all splicing junctions
            STAR --runThreadN 5 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn unaligned_reads.bam --readFilesCommand samtools view --readFilesType SAM PE --outSAMattrRGline ID:"$ID" SM:"$ID" PL:Illumina RG:"$ID" --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_2pass_"$DOWN"_

            # 2.2. MergeBamAlign
            java -jar "$PICARD" MergeBamAlignment \
                  ALIGNED="$ID"_2pass_"$DOWN"_Aligned.sortedByCoord.out.bam \
                  UNMAPPED=unaligned_reads.bam \
                  O=merge_alignments.bam \
                  R="$REFERENCE"

            # 2.3. Mark and sort duplicates
            java -jar "$PICARD" MarkDuplicates \
                  I=merge_alignments.bam \
                  O=duplicated_alignments.bam \
                  M="$ID"_"$DOWN"_marked_dup_metrics.txt

            # 2.4. Split N Cigar
            $GATK SplitNCigarReads -R "$REFERENCE" -I duplicated_alignments.bam -O splitN.bam

            # 2.5. Base Recalibration
            java -jar "$PICARD" AddOrReplaceReadGroups \
                 I=splitN.bam \
                 O=splitN.rg.bam \
                 RGID="$ID" \
                 RGLB="$ID" \
                 RGPL=Illumina \
                 RGPU=unit1 \
                 RGSM="$ID"

            $GATK BaseRecalibrator \
                  -I splitN.rg.bam \
                  -R $REFERENCE \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Homo_sapiens_assembly38.known_indels.vcf.gz" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Homo_sapiens_assembly38.dbsnp138.vcf" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
                  -O "$ID"_"$DOWN"_recal_data.table1

            $GATK ApplyBQSR \
                  -R "$REFERENCE" \
                  -I splitN.rg.bam \
                  --bqsr-recal-file "$ID"_"$DOWN"_recal_data.table1 \
                  -O "$DOWN"_recalibrated.bam
 
            $GATK BaseRecalibrator \
                  -I "$DOWN"_recalibrated.bam \
                  -R $REFERENCE \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Homo_sapiens_assembly38.known_indels.vcf.gz" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Homo_sapiens_assembly38.dbsnp138.vcf" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
                  --known-sites "/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/VCF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
                 -O "$ID"_"$DOWN"_recal_data.table2

            $GATK AnalyzeCovariates \
                  -before "$ID"_"$DOWN"_recal_data.table1 \
                  -after "$ID"_"$DOWN"_recal_data.table2 \
                  -plots "$ID"_"$DOWN"_AnalyzeCovariates.pdf

            # ALTERNATIVE: Haplotypecaller ==> Efficient but extremely slow.
            # 2.5 Variant Calling
            # $GATK --java-options "-Xmx4g" HaplotypeCaller \
            #      -R "$REFERENCE" \
            #      -I recalibrated.bam \
            #      -O "$ID"."$DOWN".vcf \
            #      -ERC BP_RESOLUTION

            # 2.5. Variant Calling
            bcftools mpileup --skip-indels --skip-indels -q 15 --threads 1 -Ou -f "$REFERENCE" "$DOWN"_recalibrated.bam | bcftools call -c -v -Ov -o  "$ID".GATK."$DOWN".vcf &

            # 2.6. Clean-up
            rm *_2pass_*_SJ.out.tab *_Log.progress.out *_Log.out 
            rm -r *__STARgenome *___STARtmp
            rm *fastq
            rm unaligned_reads.bam merge_alignments.bam splitN.bam splitN.rg.bam

      fi
done


cd /home/vincent.hahaut/data_storage/ORGANOIDS/VARIANTS/

done
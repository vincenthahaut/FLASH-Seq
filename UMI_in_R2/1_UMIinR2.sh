#!/bin/bash

ulimit -n 100000

# 0. PATH
# The sample sheet only contains the FASTQ prefix - one per lane
SAMPLESHEET="/home/vincent.hahaut/data_storage/ORGANOIDS/ID_DEEP.txt"
OUT="/home/vincent.hahaut/data_storage/ORGANOIDS/UMIinR2/"
IN="/home/vincent.hahaut/data_storage/ORGANOIDS/"

# 1. REFERENCES
STAR_75="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR/"
STAR_30="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR_29bp/"
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BBDUK/adapters.fa"

# 2. BINARIES
STAR="/home/vincent.hahaut/data_storage/binaries/STAR/bin/Linux_x86_64/STAR"
BBDUK="/home/vincent.hahaut/data_storage/binaries/bbmap/bbduk.sh"
TRIM="/home/vincent.hahaut/data_storage/binaries/Trimmomatic-0.39/trimmomatic-0.39.jar"

mkdir "$OUT"
cd "$OUT"

# 3. Load STAR References
STAR --genomeLoad LoadAndExit --genomeDir $STAR_75
STAR --genomeLoad LoadAndExit --genomeDir $STAR_30

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the CELL ID
	ID="${myArray[0]}"
	echo "$ID"
	mkdir "$ID"
	cd "$ID"

	# 1. Get the FASTQ files
	FASTQ_R1="$IN"/out/"$ID"_*_R1_001.fastq.gz
	FASTQ_R2="$IN"/out/"$ID"_*_R2_001.fastq.gz
	cat $FASTQ_R1 > sample.R1.fastq.gz
	cat $FASTQ_R2 > sample.R2.fastq.gz
	echo $FASTQ_R1
	echo $FASTQ_R2

	# 2. Get umi in R1 or R2
	if echo "$ID" | grep -qE "382_FS|316"
    then
    echo "TSO_FS_CTAAC"
	umi_tools extract --bc-pattern="^(?P<discard_1>AAGCAGTGGTATCAACGCAGAGT|AGCAGTGGTATCAACGCAGAGT|GCAGTGGTATCAACGCAGAGT|CAGTGGTATCAACGCAGAGT|AGTGGTATCAACGCAGAGT|GTGGTATCAACGCAGAGT|TGGTATCAACGCAGAGT|GGTATCAACGCAGAGT|GTATCAACGCAGAGT|ATCAACGCAGAGT|CAACGCAGAGT|AACGCAGAGT|ACGCAGAGT|GCAGAGT|CAGAGT|AGAGT|GAGT|AGT|GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=umi.UMIinR1.R1.fq --read2-in=sample.R2.fastq.gz --read2-out=umi.UMIinR1.R2.fq  --filtered-out=internal.R1.fq --filtered-out2=internal.R2.fq --extract-method=regex
 	umi_tools extract --bc-pattern="^(?P<discard_1>GAGT|AGT|GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,4})" --stdin=sample.R2.fastq.gz --stdout=umi.UMIinR2.R2.fq --read2-in=sample.R1.fastq.gz --read2-out=umi.UMIinR2.R1.fq  --filtered-out=internal.R1.fq --filtered-out2=internal.R2.fq --extract-method=regex
    elif echo "$ID" | grep -qE "385_SS3"
	then
	echo "TSO_SS3"
    umi_tools extract --bc-pattern="^(?P<discard_1>ATTGCGCAATG|TTGCGCAATG|TGCGCAATG|GCGCAATG)(?P<umi_1>.{8})GG(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=umi.UMIinR1.R1.fq --read2-in=sample.R2.fastq.gz --read2-out=umi.UMIinR1.R2.fq  --filtered-out=internal.R1.fq --filtered-out2=internal.R2.fq --extract-method=regex
    umi_tools extract --bc-pattern="^(?P<discard_1>ATTGCGCAATG|TTGCGCAATG|TGCGCAATG|GCGCAATG)(?P<umi_1>.{8})GG(?P<discard_2>G{0,4})" --stdin=sample.R2.fastq.gz --stdout=umi.UMIinR2.R2.fq --read2-in=sample.R1.fastq.gz --read2-out=umi.UMIinR2.R1.fq  --filtered-out=internal.R1.fq --filtered-out2=internal.R2.fq --extract-method=regex
    fi

	# 3. Trim reads - easier to choose the right sjdbOverhang
	java -jar $TRIM PE -threads 2 -phred33 umi.UMIinR1.R1.fq umi.UMIinR1.R2.fq umi.UMIinR1.R1.trim.fq umi.UMIinR1.R1.unpaired.fq umi.UMIinR1.R2.trim.fq umi.UMIinR1.R2.unpaired.fq CROP:75 MINLEN:25
	java -jar $TRIM PE -threads 2 -phred33 umi.UMIinR2.R1.fq umi.UMIinR2.R2.fq umi.UMIinR2.R1.trim.fq umi.UMIinR2.R1.unpaired.fq umi.UMIinR2.R2.trim.fq umi.UMIinR2.R2.unpaired.fq CROP:75 MINLEN:25

	# Only used for gene body coverage evaluation
	# Allows for a fair comparison between reads (or would be ~70x50 vs ~100x30)
	java -jar $TRIM PE -threads 2 -phred33 umi.UMIinR1.R1.fq umi.UMIinR1.R2.fq umi.UMIinR1.R1.trim30.fq umi.UMIinR1.R1.unpaired.fq umi.UMIinR1.R2.trim30.fq umi.UMIinR1.R2.unpaired.fq CROP:30 MINLEN:25
	java -jar $TRIM PE -threads 2 -phred33 umi.UMIinR2.R1.fq umi.UMIinR2.R2.fq umi.UMIinR2.R1.trim30.fq umi.UMIinR2.R1.unpaired.fq umi.UMIinR2.R2.trim30.fq umi.UMIinR2.R2.unpaired.fq CROP:30 MINLEN:25

	# 4. Zip files
	mkdir FASTQ
	gzip -c umi.UMIinR1.R1.trim.fq > FASTQ/umi.UMIinR1.R1.fq.gz
	gzip -c umi.UMIinR1.R2.trim.fq > FASTQ/umi.UMIinR1.R2.fq.gz
	gzip -c umi.UMIinR2.R1.trim.fq > FASTQ/umi.UMIinR2.R1.fq.gz
	gzip -c umi.UMIinR2.R2.trim.fq > FASTQ/umi.UMIinR2.R2.fq.gz
	#
	gzip -c umi.UMIinR1.R1.trim30.fq > FASTQ/umi.UMIinR1.R1.trim30.fq.gz
	gzip -c umi.UMIinR1.R2.trim30.fq > FASTQ/umi.UMIinR1.R2.trim30.fq.gz
	gzip -c umi.UMIinR2.R1.trim30.fq > FASTQ/umi.UMIinR2.R1.trim30.fq.gz
	gzip -c umi.UMIinR2.R2.trim30.fq > FASTQ/umi.UMIinR2.R2.trim30.fq.gz

	rm *fq

	# 5. R1-UMI: Align the Data
	mkdir STAR FEATURECOUNTS ReSQC
	if [ `zcat FASTQ/umi.UMIinR1.R1.fq.gz | grep -c ^@`  -gt 30000 ]
	then

	# 5.1. Downsampling to 30K umi reads
	seqtk sample -s42 FASTQ/umi.UMIinR1.R1.fq.gz 30000 > FASTQ/umi.30000.R1.fastq
	seqtk sample -s42 FASTQ/umi.UMIinR1.R2.fq.gz 30000 > FASTQ/umi.30000.R2.fastq
	# 5.2. MAP
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_75" --readFilesIn FASTQ/umi.30000.R1.fastq FASTQ/umi.30000.R2.fastq --seedSearchStartLmax 25 --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_UMIinR1_downsampled_
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.out.bam > STAR/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.filtered.bam
	# 5.3. READ DISTIRBUTION
	read_distribution.py -i STAR/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_UMIinR1_readDistribution.txt &
	# 5.4. ASSIGN FEATURES
	featureCounts -T 1 -p -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_UMIinR1_30000.feature.txt STAR/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.filtered.bam
	samtools view -h FEATURECOUNTS/"$ID"_UMIinR1_downsampled_Aligned.sortedByCoord.filtered.bam.featureCounts.bam | awk -F'\t' '$1~/_.{8}$|@SQ|@HD/{print $0}' | samtools view -Sb - > "$ID"_assignedUMIinR1_30000.bam 
	samtools sort -@ 10 "$ID"_assignedUMIinR1_30000.bam  > "$ID"_assignedUMIinR1_30000.sorted.bam
	samtools index "$ID"_assignedUMIinR1_30000.sorted.bam
	# 5.5. COUNT UMI
	umi_tools count --per-gene --paired --gene-tag=XT --chimeric-pairs=discard --unpaired-reads=discard --assigned-status-tag=XS -I "$ID"_assignedUMIinR1_30000.sorted.bam -S FEATURECOUNTS/"$ID".umi.UMIinR1.30000.tsv.gz
	fi

	# 6. R2-UMI: Align the Data
	if [ `zcat FASTQ/umi.UMIinR2.R1.fq.gz | grep -c ^@`  -gt 30000 ]
	then
	
	# 6.1. Downsampling to 30K umi reads
	seqtk sample -s42 FASTQ/umi.UMIinR2.R1.fq.gz 30000 > FASTQ/umi.30000.R1.fastq
	seqtk sample -s42 FASTQ/umi.UMIinR2.R2.fq.gz 30000 > FASTQ/umi.30000.R2.fastq
	# 6.2. MAP
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_75" --readFilesIn FASTQ/umi.30000.R2.fastq FASTQ/umi.30000.R1.fastq --seedSearchStartLmax 25 --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_UMIinR2_downsampled_
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.out.bam > STAR/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.filtered.bam
	# 6.3. READ DISTRIBUTION
	read_distribution.py -i STAR/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_UMIinR2_readDistribution.txt &
	# 6.4. ASSIGN FEATURES
	featureCounts -T 1 -p -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_UMIinR2_30000.feature.txt STAR/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.filtered.bam
	samtools view -h FEATURECOUNTS/"$ID"_UMIinR2_downsampled_Aligned.sortedByCoord.filtered.bam.featureCounts.bam | awk -F'\t' '$1~/_.{8}$|@SQ|@HD/{print $0}' | samtools view -Sb - > "$ID"_assignedUMIinR2_30000.bam 
	samtools sort -@ 10 "$ID"_assignedUMIinR2_30000.bam  > "$ID"_assignedUMIinR2_30000.sorted.bam
	samtools index "$ID"_assignedUMIinR2_30000.sorted.bam
	# 6.5. COUNT UMI
	umi_tools count --per-gene --paired --gene-tag=XT --chimeric-pairs=discard --unpaired-reads=discard --assigned-status-tag=XS -I "$ID"_assignedUMIinR2_30000.sorted.bam -S FEATURECOUNTS/"$ID".umi.UMIinR2.30000.tsv.gz
	fi
	
	# 7.1. Align 30x30bp for gene body coverage comparison
	# need to adapt seedSearch
	# After a few test --seedSearchStartLmax 25 works well
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_30" --readFilesIn FASTQ/umi.UMIinR1.R1.trim30.fq.gz FASTQ/umi.UMIinR1.R2.trim30.fq.gz --seedSearchStartLmax 25 --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_UMIinR1_trimmed_
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_UMIinR1_trimmed_Aligned.sortedByCoord.out.bam > STAR/"$ID"_UMIinR1_trimmed_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_UMIinR1_trimmed_Aligned.sortedByCoord.filtered.bam
	#
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_30" --readFilesIn FASTQ/umi.UMIinR2.R2.trim30.fq.gz FASTQ/umi.UMIinR2.R1.trim30.fq.gz --seedSearchStartLmax 25 --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_UMIinR2_trimmed_
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_UMIinR2_trimmed_Aligned.sortedByCoord.out.bam > STAR/"$ID"_UMIinR2_trimmed_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_UMIinR2_trimmed_Aligned.sortedByCoord.filtered.bam
	
	# 8. ReSQC
	mkdir ReSQC
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_UMIinR1_trimmed_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_UMIinR1_genebody &
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_UMIinR2_trimmed_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_UMIinR2_genebody &



	cd "$OUT"
done

# Unload STAR indexes
STAR --genomeLoad Remove --genomeDir "$STAR_75"
STAR --genomeLoad Remove --genomeDir "$STAR_30"

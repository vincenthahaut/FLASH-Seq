#!/bin/bash

ulimit -n 100000

# 0. PATH
# The sample sheet only contains the FASTQ prefix - one per lane
SAMPLESHEET="/home/vincent.hahaut/data_storage/ORGANOIDS/316_ae"
OUT="/home/vincent.hahaut/data_storage/ORGANOIDS/STAR_DEEP/"
IN="/home/vincent.hahaut/data_storage/ORGANOIDS/"

# 1. REFERENCES
STAR_100="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR_100bp/"
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
GTF_EXONINTRON="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.exonsANDintrons.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BBDUK/adapters.fa"
SALMON_ALL="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/salmon_precomputed/salmon_sa_index/"
FASTA="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/DNA/GRCh38.primary_assembly.genome.fa"

# 2. BINARIES
STAR="/home/vincent.hahaut/data_storage/binaries/STAR/bin/Linux_x86_64/STAR"
BBDUK="/home/vincent.hahaut/data_storage/binaries/bbmap/bbduk.sh"
BBSPLIT="/home/vincent.hahaut/data_storage/binaries/bbmap/bbsplit.sh"
BBMAP_filter="/home/vincent.hahaut/data_storage/binaries/bbmap/filterbyname.sh"
bbmapsaturation="/home/vincent.hahaut/data_storage/binaries/bbmap/bbcountunique.sh"
TRIM="/home/vincent.hahaut/data_storage/binaries/Trimmomatic-0.39/trimmomatic-0.39.jar"
salmon="/home/vincent.hahaut/data_storage/binaries/salmon-1.5.2_linux_x86_64/bin/salmon"
filerInvasion="/home/vincent.hahaut/Desktop/FLASH-Seq/R/UMI/strand_invasion/filterInvasionEventsfromBAM.R"

mkdir "$OUT"
cd "$OUT"

# 3. Load STAR References
STAR --genomeLoad LoadAndExit --genomeDir $STAR_100

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the CELL ID
	ID="${myArray[0]}"
	echo "$ID"
	mkdir "$ID"
	cd "$ID"

	# 1. Get the FASTQ R1
	FASTQ_R1="$IN"/out/"$ID"_*_R1_001.fastq.gz
	FASTQ_R2="$IN"/out/"$ID"_*_R2_001.fastq.gz
	cat $FASTQ_R1 > sample.R1.fastq.gz
	cat $FASTQ_R2 > sample.R2.fastq.gz
	echo $FASTQ_R1
	echo $FASTQ_R2

	# 2. Get UMI
	# Depending on the cell ID - extract the right UMI - Spacer sequence
	if echo "$ID" | grep -qE "FS_315|315_FS|316_FS|382_FS|317_FS|318_FS"
	then
	echo "TSO_FS_CTAAC"
	umi_tools extract --bc-pattern="^(?P<discard_1>AAGCAGTGGTATCAACGCAGAGT|AGCAGTGGTATCAACGCAGAGT|GCAGTGGTATCAACGCAGAGT|CAGTGGTATCAACGCAGAGT|AGTGGTATCAACGCAGAGT|GTGGTATCAACGCAGAGT|TGGTATCAACGCAGAGT|GGTATCAACGCAGAGT|GTATCAACGCAGAGT|ATCAACGCAGAGT|CAACGCAGAGT|AACGCAGAGT|ACGCAGAGT|GCAGAGT|CAGAGT|AGAGT|GAGT|AGT|GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=umi.UMIinR1.R1.fq --read2-in=sample.R2.fastq.gz --read2-out=umi.UMIinR1.R2.fq --extract-method=regex
	umi_tools extract --bc-pattern="^(?P<discard_1>GAGT|AGT|GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,4})" --stdin=sample.R2.fastq.gz --stdout=umi.UMIinR2.R2.fq --read2-in=sample.R1.fastq.gz --read2-out=umi.UMIinR2.R1.fq --extract-method=regex
	# In very rare cases (<0.0001%) can get the UMI in both R1 and R2
	# Find reads with a UMI in both R1 and R2
	cat umi.UMIinR1.R1.fq | uniq | awk 'NR%4==1{print}' | sed 's/\@//g' > names.R1umi.txt
	cat umi.UMIinR2.R1.fq | uniq | awk 'NR%4==1{print}' | sed 's/\@//g' > names.R2umi.txt
	# Remove the UMI info from their name
	cat names.R1umi.txt | sed 's/_........ .*$//g' > names.R1umi.cleaned
	cat names.R2umi.txt | sed 's/_........ .*$//g' > names.R2umi.cleaned
	comm -12 <(sort names.R1umi.cleaned) <(sort names.R2umi.cleaned) > R1R2.toFilterOut
	echo "===> Number of R1-R2 with both a UMI: $(wc -l R1R2.toFilterOut) <==="
	echo "===> Number of R1 UMI before cleanup: $(wc -l names.R1umi.txt) <==="
	echo "===> Number of R2 UMI before cleanup: $(wc -l names.R2umi.txt) <==="
	# Get their full ID
	grep -f R1R2.toFilterOut names.R1umi.txt > R1.toFilterOut
	grep -f R1R2.toFilterOut names.R2umi.txt > R2.toFilterOut
	# Remove them from the data to get clean UMI files
	$BBMAP_filter -Xmx6g in=umi.UMIinR1.R1.fq in2=umi.UMIinR1.R2.fq out=umi.UMIinR1.R1.tmp out2=umi.UMIinR1.R2.tmp names=R1.toFilterOut include=f overwrite=t
	$BBMAP_filter -Xmx6g in=umi.UMIinR2.R1.fq in2=umi.UMIinR2.R2.fq out=umi.UMIinR2.R1.tmp out2=umi.UMIinR2.R2.tmp names=R2.toFilterOut include=f overwrite=t
	mv umi.UMIinR1.R1.tmp umi.UMIinR1.R1.fq
	mv umi.UMIinR2.R1.tmp umi.UMIinR2.R1.fq
	mv umi.UMIinR1.R2.tmp umi.UMIinR1.R2.fq
	mv umi.UMIinR2.R2.tmp umi.UMIinR2.R2.fq
	echo "===> Number of R1 UMI after cleanup: $(grep -c \@ umi.UMIinR1.R1.fq) <==="
	echo "===> Number of R2 UMI after cleanup: $(grep -c \@ umi.UMIinR2.R1.fq) <==="
	# Get Internal reads only by excluding the UMI reads
	cat umi.UMIinR1.R1.fq | uniq | awk 'NR%4==1{print}' | sed 's/\@//g' | sed 's/_........//g' > names.R1umi.txt
	cat umi.UMIinR2.R1.fq | uniq | awk 'NR%4==1{print}' | sed 's/\@//g' | sed 's/_........//g' > names.R2umi.txt
	sed -i 's/_........//g' R1.toFilterOut
	sed -i 's/_........//g' R2.toFilterOut
	cat names.R1umi.txt names.R2umi.txt R1.toFilterOut R2.toFilterOut > names.umi.txt
	$BBMAP_filter -Xmx6g in=sample.R1.fastq.gz in2=sample.R2.fastq.gz out=internal.R1.fq out2=internal.R2.fq names=names.umi.txt include=f overwrite=t
	rm R1.toFilterOut R2.toFilterOut toFilterOut.txt names.umi.txt
	echo "===> Number of Internal Reads after cleanup: $(grep -c \@ internal.R1.fq) <==="
	# Combine UMI reads.
	# R2-UMI are in opposite direction compared to the gene. To reconciliate them use R2 R2-UMI as R1.
	cat umi.UMIinR1.R1.fq umi.UMIinR2.R2.fq > umi.R1.fq
	cat umi.UMIinR1.R2.fq umi.UMIinR2.R1.fq > umi.R2.fq
	# Clean-up
	rm names.* umi.UMIinR*.R*.fq
	# Not used:
	# elif echo "$ID" | grep -qE "SS3_333|385_SS3"
	# then
	# echo "TSO_SS3"
	# umi_tools extract --bc-pattern="^(?P<discard_1>ATTGCGCAATG){s<=1}(?P<umi_1>.{8})GG(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=umi.R1.fq --read2-in=sample.R2.fastq.gz --read2-out=umi.R2.fq  --filtered-out=internal.R1.fq --filtered-out2=internal.R2.fq --extract-method=regex
	fi

	mkdir FASTQ 
	cat umi.R1.fq internal.R1.fq > FASTQ/allreads.R1.fq
	cat umi.R2.fq internal.R2.fq  > FASTQ/allreads.R2.fq
	mv umi.R*.fq FASTQ/
	mv internal.R*.fq FASTQ/

    # 3. Trim Reads
  	# t = threads
	# ktrim = right and left
	# rcomp=f only look for the forward sequence not reverse
	# k = min kmer length trim
	# hdist = hamming distance for mismatches
	# mink = min kmer length trim - at edges of fragments
	# hdist2 = for kmer < k, use this hamming distance
	# minlength = minimal read length or discard
	# tbo = trim on where paired reads overlap 
	# Remove F/R sequencing adapters
    $BBDUK -Xmx6g in1=FASTQ/allreads.R1.fq in2=FASTQ/allreads.R2.fq out1=FASTQ/allreads.R1.2.fq out2=FASTQ/allreads.R2.2.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    $BBDUK -Xmx6g in1=FASTQ/allreads.R1.2.fq in2=FASTQ/allreads.R2.2.fq out1=FASTQ/allreads.R1.trim.fq out2=FASTQ/allreads.R2.trim.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo

    $BBDUK -Xmx6g in1=FASTQ/umi.R1.fq in2=FASTQ/umi.R2.fq out1=FASTQ/umi.R1.2.fq out2=FASTQ/umi.R2.2.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    $BBDUK -Xmx6g in1=FASTQ/umi.R1.2.fq in2=FASTQ/umi.R2.2.fq out1=FASTQ/umi.R1.trim.fq out2=FASTQ/umi.R2.trim.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo

    $BBDUK -Xmx6g in1=FASTQ/internal.R1.fq in2=FASTQ/internal.R2.fq out1=FASTQ/internal.R1.2.fq out2=FASTQ/internal.R2.2.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    $BBDUK -Xmx6g in1=FASTQ/internal.R1.2.fq in2=FASTQ/internal.R2.2.fq out1=FASTQ/internal.R1.trim.fq out2=FASTQ/internal.R2.trim.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo

	rm FASTQ/allreads.R1.2.fq FASTQ/internal.R1.2.fq FASTQ/umi.R1.2.fq FASTQ/allreads.R2.2.fq FASTQ/internal.R2.2.fq FASTQ/umi.R2.2.fq

	# 4. Discard (and trim if needed) reads
	# After testing a few different mapping strategies to accomodate with the varying read length, it appears that a long sjdbOverhang (100) and a shorter seedSearchStartLmax (30) works the best for UMI reads in terms of % mapping and canonical splice junctions
	java -jar $TRIM PE -threads 2 -phred33 FASTQ/umi.R1.trim.fq FASTQ/umi.R2.trim.fq FASTQ/umi.R1.trim2.R1.fq FASTQ/umi.R1.trim2.unpaired.fq FASTQ/umi.R2.trim2.R1.fq FASTQ/umi.R2.trim2.unpaired.fq CROP:100 MINLEN:29
	mv FASTQ/umi.R1.trim2.R1.fq FASTQ/umi.R1.trim.fq
	mv FASTQ/umi.R2.trim2.R1.fq FASTQ/umi.R2.trim.fq

	# 5. Gather the FASTQ 
	gzip -c FASTQ/allreads.R1.trim.fq > FASTQ/allreads.R1.fq.gz
	gzip -c FASTQ/umi.R1.trim.fq > FASTQ/umi.R1.fq.gz
	gzip -c FASTQ/internal.R1.trim.fq > FASTQ/internal.R1.fq.gz
	rm FASTQ/allreads.R1.trim.fq FASTQ/umi.R1.trim.fq FASTQ/internal.R1.trim.fq
	#
	gzip -c FASTQ/allreads.R2.trim.fq > FASTQ/allreads.R2.fq.gz
	gzip -c FASTQ/umi.R2.trim.fq > FASTQ/umi.R2.fq.gz
	gzip -c FASTQ/internal.R2.trim.fq > FASTQ/internal.R2.fq.gz
	rm FASTQ/allreads.R2.trim.fq FASTQ/umi.R2.trim.fq FASTQ/internal.R2.trim.fq
	
	# 6. Redefine variables
	MYFASTQ_R1="$OUT"/"$ID"/FASTQ/allreads.R1.fq.gz
	MYFASTQ_UMI_R1="$OUT"/"$ID"/FASTQ/umi.R1.fq.gz
	MYFASTQ_INTERNAL_R1="$OUT"/"$ID"/FASTQ/internal.R1.fq.gz
	MYFASTQ_R2="$OUT"/"$ID"/FASTQ/allreads.R2.fq.gz
	MYFASTQ_UMI_R2="$OUT"/"$ID"/FASTQ/umi.R2.fq.gz
	MYFASTQ_INTERNAL_R2="$OUT"/"$ID"/FASTQ/internal.R2.fq.gz

	if [ `zcat "$MYFASTQ_R1" | grep -c ^@`  -gt "25000" ]
	then

	# 7. Align the Data
	mkdir STAR

	# 7.1. Untrimmed reads - For debugging purposes
	# OUTPUT EVERYTHING --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outReadsUnmapped Fastx
	# $STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_100" --readFilesIn sample.R1.fastq.gz sample.R2.fastq.gz --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_ALLnoTrim_
	# samtools view -@ 30 -Sb -F 260 -q 5 STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.out.bam > STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam
	# samtools index STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam

	# 7.2. Trimmed reads
	for READS in "$MYFASTQ_R1" "$MYFASTQ_UMI_R1" "$MYFASTQ_INTERNAL_R1"
	do

	# 7.2.1. Set Variables
	if echo "$READS" | grep -q "allreads"
	then
	TYPE="ALL"
	R1="$MYFASTQ_R1"
	R2="$MYFASTQ_R2"
	elif echo "$READS" | grep -q "umi"
	then
	TYPE="UMI"
	R1="$MYFASTQ_UMI_R1"
	R2="$MYFASTQ_UMI_R2"
	elif echo "$READS" | grep -q "internal"
	then
	TYPE="INTERNAL"
	R1="$MYFASTQ_INTERNAL_R1"
	R2="$MYFASTQ_INTERNAL_R2"
	fi

	echo "$READS"
	echo $TYPE	

	# 7.2.2. Map and count reads: All reads or internal reads
	mkdir FEATURECOUNTS
	if [ "$TYPE" != "UMI" ]
	then
	# 7.2.2.1. Alignment
	$STAR --runThreadN 10 --limitBAMsortRAM 10000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_100" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_"$TYPE"_
	# 7.2.2.2. Filter out unmapped and secundary alignments
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	# 7.2.2.3. Count
	featureCounts -T 1 -p -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.GENE.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	featureCounts -T 1 -p -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.EXON.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	# 7.2.2.4. Transcript count
	if [ "$TYPE" == "ALL" ]
	then
	$salmon quant -i $SALMON_ALL -p 4 -l IU --minAssignedFrags 1 --dumpEqWeights -1 $R1 -2 $R2 --validateMappings -o STAR/"$ID"_salmon_ALL_TRANSCRIPTS
	fi
	fi

	# 7.2.3. Map and count reads: UMI reads
	if [ "$TYPE" == "UMI" ]
	then
	# 7.2.3.1. Alignment
	$STAR --runThreadN 10 --limitBAMsortRAM 10000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_100" --readFilesIn "$R1" "$R2" --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --seedSearchStartLmax 30 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_"$TYPE"_
	# 7.2.3.2. Filter out unmapped and secundary alignments
	samtools view -@ 10 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	# 7.2.3.3. Filter out invading UMI
	# arg1: path/to/bam
  	# arg2: sample ID
  	# arg3: method ('FLASH-Seq' or 'zUMIs')
  	# arg4: number of mistmaches (recommended: 1)
  	# arg5: path/to/genome.fa - genome used for the mapping
  	# arg6: minimum number of total reads to process the sample (recommended: 10000)
  	# arg7: minimum mapq to look for a match between UMI and adjacent sequence (recommended: 5)
  	# arg8: path/to/output.bam"
  	Rscript $filerInvasion STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam "$ID"_"$TYPE" 'FLASH-Seq' 1 $FASTA 1000 5 STAR/
  	mv STAR/"$ID"_UMI_filteredInvasion.bam STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	# 7.2.3.4. Count
	featureCounts -T 1 -p -t gene -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_Aligned.onGENE.txt STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam &	
	featureCounts -T 1 -p -t all -g feature -s 0 --fracOverlap 0.25 -a "$GTF_EXONINTRON" -R BAM -o FEATURECOUNTS/"$ID"_unstranded._UMI_featureCounts_ExonIntrons.txt STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam
	mv FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.bam FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.exonIntron.bam
	# 7.2.3.5. Assign reads to a feature
	featureCounts -T 1 -p -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_Aligned.onEXON.txt STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam
	samtools sort  -@ 10 FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.bam -o FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.sorted.bam
	samtools index FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.sorted.bam
	rm FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.bam
	# 7.2.3.6. Count UMI reads
	umi_tools dedup --per-gene --paired --gene-tag=XT --chimeric-pairs=discard --unpaired-reads=discard --assigned-status-tag=XS -I FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.sorted.bam -S FEATURECOUNTS/"$ID".umi.dedup.bam --output-stats=FEATURECOUNTS/"$ID".umi.dedup.
	umi_tools count --per-gene --paired --gene-tag=XT --chimeric-pairs=discard --unpaired-reads=discard --assigned-status-tag=XS -I FEATURECOUNTS/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam.featureCounts.sorted.bam -S FEATURECOUNTS/"$ID".umi.counts.tsv.gz
	# 7.2.3.7. Transcript count
	$salmon quant -i $SALMON_ALL -p 4 -l ISF --minAssignedFrags 1 --dumpEqWeights -1 $R1 -2 $R2 --validateMappings -o STAR/"$ID"_salmon_UMI_TRANSCRIPTS
	fi

	# 8. Downsampling
	mkdir DOWNSAMPLE
	cd DOWNSAMPLE
	for DOWN in 5000 10000 20000 40000 50000 75000 100000 125000 150000 200000 300000 400000 500000
	do
	if [ `zcat "$R1" | grep -c ^@`  -gt "$DOWN" ]
	then
	mkdir DOWN_"$DOWN"
	cd DOWN_"$DOWN"
	mkdir FEATURECOUNTS
	# 8.1. Downsample reads
	seqtk sample -s42 "$R1" "$DOWN" > "$TYPE"."$DOWN".R1.fastq
	seqtk sample -s42 "$R2" "$DOWN" > "$TYPE"."$DOWN".R2.fastq

	#### 8.2. UMI ####
	if [ "$TYPE" == "UMI" ]
	then
	# 8.2.1. Map and count reads: UMI reads
	$STAR --runThreadN 10 --limitBAMsortRAM 10000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_100" --readFilesIn "$TYPE"."$DOWN".R1.fastq "$TYPE"."$DOWN".R2.fastq --readFilesCommand cat --seedSearchStartLmax 30 --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_"$TYPE"_"$DOWN"_
	samtools view -@ 10 -Sb -F 260 -q 5 "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam > "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	rm "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam
	# 8.2.2. Filter out invading UMI
  	Rscript $filerInvasion "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam "$ID"_"$TYPE"_"$DOWN" 'FLASH-Seq' 1 $FASTA 1000 5 ./
  	mv "$ID"_"$TYPE"_"$DOWN"_filteredInvasion.bam "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	# 8.2.3. Assign reads to a feature
	featureCounts -T 5 -p -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools view -h FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam | awk -F'\t' '$1~/_.{8}$|@SQ|@HD/{print $0}' | samtools view -Sb - > "$ID"_assignedUMI_"$DOWN".bam 
	samtools sort -@ 10 "$ID"_assignedUMI_"$DOWN".bam  > "$ID"_assignedUMI_"$DOWN".sorted.bam
	rm FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam
	samtools index "$ID"_assignedUMI_"$DOWN".sorted.bam
	# 8.2.3. Count UMI
	umi_tools count --per-gene --paired --gene-tag=XT --chimeric-pairs=discard --unpaired-reads=discard --assigned-status-tag=XS -I "$ID"_assignedUMI_"$DOWN".sorted.bam -S FEATURECOUNTS/"$ID".umi.counts."$DOWN".tsv.gz
	# 8.2.4. Transcript count
	$salmon quant -i $SALMON_ALL -p 4 -l ISF --minAssignedFrags 1 --dumpEqWeights -1 "$TYPE"."$DOWN".R1.fastq -2 "$TYPE"."$DOWN".R2.fastq --validateMappings -o STAR/"$ID"_salmon_UMI_"$DOWN"_TRANSCRIPTS &

	#### 8.3. INTERNAL ####
	elif [ "$TYPE" != "UMI" ]
	then
	# 8.3.1. Map and count reads: INTERNAL reads
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_100" --readFilesIn "$TYPE"."$DOWN".R1.fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_"$TYPE"_"$DOWN"_
	samtools view -@ 10 -Sb -F 260 -q 5 "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam > "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	rm "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam
	# 8.3.2. Assign reads to a feature
	featureCounts -T 1 -p -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.GENE.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	featureCounts -T 1 -p -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.EXON.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	# 8.3.3. Transcript count
	$salmon quant -i $SALMON_ALL -p 4 -l IU --minAssignedFrags 1 --dumpEqWeights -1 "$TYPE"."$DOWN".R1.fastq -2 "$TYPE"."$DOWN".R2.fastq --validateMappings -o STAR/"$ID"_salmon_INTERNAL_"$DOWN"_TRANSCRIPTS &
	fi

	cd "$OUT"/"$ID"/DOWNSAMPLE/

	fi

	done

	cd "$OUT"/"$ID"/

	done

	# 9. ReSQC
	cd "$OUT"/"$ID"/
	mkdir ReSQC
	# 9.1. GeneBodyCoverage
	# Heavy slow down
	# geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_INTERNAL_geneBody.all &
	# geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_UMI_geneBody.all &
	# 9.2. ReSQC Statistics
	read_distribution.py -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_INTERNAL_readDistribution.txt &
	read_distribution.py -i STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_UMI_readDistribution.txt &

    fi

	cd "$OUT"
done

# Unload STAR indexes
STAR --genomeLoad Remove --genomeDir "$STAR_100"

#!/bin/bash

ulimit -n 100000

# 0. PATH
# The sample sheet only contains the FASTQ prefix - one per lane
SAMPLESHEET="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/xad"
OUT="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/"
IN="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/"

# 1. REFERENCES
STAR_75="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR/" # Internal Reads
STAR_50="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR_50bp/" # UMI Reads
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
GTF_EXONINTRON="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.exonsANDintrons.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BBDUK/adapters.fa"
SALMON_ALL="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/salmon_precomputed/salmon_sa_index/"

# 2. BINARIES
STAR="/home/vincent.hahaut/data_storage/binaries/STAR/bin/Linux_x86_64/STAR"
BBDUK="/home/vincent.hahaut/data_storage/binaries/bbmap/bbduk.sh"
BBMAP_filter="/home/vincent.hahaut/data_storage/binaries/bbmap/filterbyname.sh"
TRIM="/home/vincent.hahaut/data_storage/binaries/Trimmomatic-0.39/trimmomatic-0.39.jar"
salmon=/home/vincent.hahaut/data_storage/binaries/salmon-1.5.2_linux_x86_64/bin/salmon

mkdir "$OUT"
cd "$OUT"

# 3. Load STAR References
STAR --genomeLoad LoadAndExit --genomeDir $STAR_75
STAR --genomeLoad LoadAndExit --genomeDir $STAR_50

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the CELL ID
	ID="${myArray[0]}"
	echo "$ID"	
	mkdir "$ID"
	cd "$ID"

	# 1. Get the FASTQ R1
	#FASTQ_R1=$(echo $IN | tr "," "\n" | awk -v ID=$ID '{print $1"/out/"ID"*_R1_*fastq.gz"}')
	FASTQ_R1="$IN"/out/"$ID"*_R1_*.fastq.gz
	cat $FASTQ_R1 > sample.R1.fastq.gz

	# 2. Get UMI
	# Depending on the cell ID - extract the right UMI - Spacer sequence
	mkdir FASTQ 

    if echo "$ID" | grep -q "Spacer1"
    then
    echo "TSO-Spacer1-CAGCA"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CAGCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -qE "SS3fwd"
	then
	echo "SS3-TSO"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*ATTGCGCAATG){s<=1}(?P<umi_1>.{8})(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "Spacer2"
    then
    echo "TSO-Spacer2-ATAAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>ATAACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
    #
    elif echo "$ID" | grep -q "TSO_SS3_CAGCA"
    then
    echo "TSO_SS3_CAGCA"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*TG)(?P<umi_1>.{8})(?P<discard_2>CAGCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
    #
    elif echo "$ID" | grep -qE "TSO_SS3_CGTAC|311_SS3"
    then
    echo "TSO_SS3_CAGCA"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*TG)(?P<umi_1>.{8})(?P<discard_2>CGTACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
    #
    elif echo "$ID" | grep -q "TSO_SS3_CTGAC"
    then
    echo "TSO_SS3_CTGAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*TG)(?P<umi_1>.{8})(?P<discard_2>CTGACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
    #
	elif echo "$ID" | grep -qE "CTAAC|311_CTAAC"
    then
    echo "TSO-CTAAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "TSO_FS_CGTAC"
    then
    echo "TSO-CGTAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CGTACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "ATGAC_"
    then
    echo "TSO-ATGAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>ATGACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "CTGAC_"
    then
    echo "TSO-CTGAC|TSO_SS3_CTGAC"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CTGACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "AAGCA_"
    then
    echo "TSO-AAGCA"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>AAGCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -q "CATCA_"
    then
    echo "TSO-CATCA"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CATCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	#
	elif echo "$ID" | grep -qE "FStso|Fstso"
	then
	# Could likely reduce the motif length but would start to be more dangerous. 
	echo "FS-TSO"
   	umi_tools extract --bc-pattern="^(?P<discard_1>.+CAACGCAGAGT){s<=1}(?P<umi_1>.{8})(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -qE "5_SS3|6_SS3|277_TSO_SS3|304_TSO_SS3"
	then
	echo "SS3-TSO"
    umi_tools extract --bc-pattern="^(?P<discard_1>.*ATTGCGCAATG){s<=1}(?P<umi_1>.{8})(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	fi

	# 3. UMI-tools has a filtered-out option that could get the internal reads out but it does not work here. Seems buggy with single-end reads.
	# Extract Internal reads
	awk 'NR%4==1{print}' FASTQ/umi.R1.fq | sed 's/_........//g' | sed 's/^@//g' > names.txt
	$BBMAP_filter in=sample.R1.fastq.gz out=FASTQ/internal.R1.fq names=names.txt include=f overwrite=t
	rm names.txt

	cat FASTQ/umi.R1.fq FASTQ/internal.R1.fq > FASTQ/allreads.R1.fq

	# 4. Remove leftover of adapters
    # Remove F/R sequencing adapters
  	# t = threads
	# ktrim = right and left
	# rcomp=f only look for the forward sequence not reverse
	# k = min kmer length trim
	# hdist = hamming distance for mismatches
	# mink = min kmer length trim - at edges of fragments
	# hdist2 = for kmer < k, use this hamming distance
	# minlength = minimal read length or discard
	# tbo = trim on where paired reads overlap
	$BBDUK -Xmx48g in=FASTQ/allreads.R1.fq out=FASTQ/allreads.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 minlength=30 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/allreads.R1.2.fq out=FASTQ/allreads.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 minlength=30 mink=10 hdist=1
  	
  	$BBDUK -Xmx48g in=FASTQ/umi.R1.fq out=FASTQ/umi.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=18 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/umi.R1.2.fq out=FASTQ/umi.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=18 mink=10 hdist=1

  	$BBDUK -Xmx48g in=FASTQ/internal.R1.fq out=FASTQ/internal.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/internal.R1.2.fq out=FASTQ/internal.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=10 hdist=1

	rm FASTQ/allreads.R1.2.fq FASTQ/internal.R1.2.fq FASTQ/umi.R1.2.fq

	# 5. Different read length were used. Trim everything to the same length
	java -jar $TRIM SE -threads 10 -phred33 FASTQ/internal.R1.trim.fq FASTQ/internal.R1.trim2.fq CROP:75 MINLEN:60
	java -jar $TRIM SE -threads 10 -phred33 FASTQ/umi.R1.trim.fq FASTQ/umi.R1.trim2.R1.fq CROP:50 MINLEN:30
	mv FASTQ/internal.R1.trim2.fq FASTQ/internal.R1.trim.fq
	mv FASTQ/umi.R1.trim2.R1.fq FASTQ/umi.R1.trim.fq

	# 6. Gather the FASTQ and redefine the variables
	gzip -c FASTQ/allreads.R1.trim.fq > FASTQ/allreads.R1.fq.gz
	gzip -c FASTQ/umi.R1.trim.fq > FASTQ/umi.R1.fq.gz
	gzip -c FASTQ/internal.R1.trim.fq > FASTQ/internal.R1.fq.gz
	#
	MYFASTQ_R1="$OUT"/"$ID"/FASTQ/allreads.R1.fq.gz
	MYFASTQ_UMI_R1="$OUT"/"$ID"/FASTQ/umi.R1.fq.gz
	MYFASTQ_INTERNAL_R1="$OUT"/"$ID"/FASTQ/internal.R1.fq.gz
	#
	rm FASTQ/allreads.R1.trim.fq FASTQ/umi.R1.trim.fq FASTQ/internal.R1.trim.fq

	# 7. Align the Data
	mkdir STAR

	# 7.1. Untrimmed reads - For debugging purposes
	# OUTPUT EVERYTHING --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outReadsUnmapped Fastx
	# $STAR --runThreadN 5 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_75" --readFilesIn sample.R1.fastq.gz --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_ALLnoTrim_
	# samtools view -@ 5 -Sb -F 260 -q 5 STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.out.bam > STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam
	# samtools index STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam

	# 7.2. Trimmed reads
	for READS in "$MYFASTQ_R1" "$MYFASTQ_UMI_R1" "$MYFASTQ_INTERNAL_R1"
	do

	# 7.2.1. Set Variables
	if echo "$READS" | grep -q "allreads"
	then
	TYPE="ALL"
	elif echo "$READS" | grep -q "umi"
	then
	TYPE="UMI"
	elif echo "$READS" | grep -q "internal"
	then
	TYPE="INTERNAL"
	fi

	echo "$READS"
	echo $TYPE	

	# 7.2.2. Map and count reads: All reads or internal reads
	mkdir FEATURECOUNTS
	if [ "$TYPE" != "UMI" ]
	then
	# 7.2.2.1. Alignment
	$STAR --runThreadN 5 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_75" --readFilesIn "$READS" --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_"$TYPE"_
	# 7.2.2.2. Filter out unmapped and secundary alignments
	samtools view -@ 5 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	# 7.2.2.3. Count
	featureCounts -T 1 -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.GENE.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	featureCounts -T 1 -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.EXON.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	# 7.2.2.4. Transcript count
	if [ "$TYPE" == "ALL" ]
	then
	$salmon quant -i $SALMON_ALL -p 5 -l U --fldMean 700 --fldSD 200 --fldMax 1100 --minAssignedFrags 1 --dumpEqWeights -r $READS --validateMappings -o STAR/"$ID"_salmon_ALL_TRANSCRIPTS
	fi
	fi

	# 7.2.3. Map and count reads: UMI reads
	if [ "$TYPE" == "UMI" ]
	then
	# 7.2.3.1. Alignment
	$STAR --runThreadN 5 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_50" --readFilesIn "$READS" --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_"$TYPE"_
	# 7.2.3.2. Filter out unmapped and secundary alignments
	samtools view -@ 5 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	# 7.2.3.3. Count
	mkdir FEATURECOUNTS/GENE FEATURECOUNTS/EXON FEATURECOUNTS/EXONINTRON
	featureCounts -T 1 -t gene -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM  --Rpath FEATURECOUNTS/GENE/ -o FEATURECOUNTS/"$ID"_GENE.STRANDED STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam &
	featureCounts -T 1 -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM  --Rpath FEATURECOUNTS/EXON/ -o FEATURECOUNTS/"$ID"_EXON.STRANDED STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam &
	featureCounts -T 1 -t all -g feature_name -s 0 --largestOverlap --fracOverlap 0.25 -a "$GTF_EXONINTRON" -R BAM --Rpath FEATURECOUNTS/EXONINTRON/ -o FEATURECOUNTS/"$ID"_UMI_featureCounts_ExonIntrons STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam & 
	# 7.2.3.4. Transcript count
	$salmon quant -i $SALMON_ALL -p 5 -l SF --fldMean 700 --fldSD 200 --fldMax 1100 --minAssignedFrags 1 --dumpEqWeights -r $READS --validateMappings -o STAR/"$ID"_salmon_UMI_TRANSCRIPTS
	cd "$OUT"
	fi

	# 8. Downsampling
	mkdir DOWNSAMPLE
	cd DOWNSAMPLE
	for DOWN in 5000 10000 20000 40000 50000 75000 100000 125000 250000
	do
	if [ `zcat "$READS" | grep -c ^@`  -gt "$DOWN" ]
	then
	mkdir DOWN_"$DOWN"
	cd DOWN_"$DOWN"
	mkdir FEATURECOUNTS
	# 8.1. Downsample reads
	seqtk sample -s42 "$READS" "$DOWN" > "$TYPE"."$DOWN".R1.fastq
	
	#### 8.2. UMI ####
	if [ "$TYPE" == "UMI" ]
	then
	# 8.2.1. Map and count reads: UMI reads
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_50" --readFilesIn "$TYPE"."$DOWN".R1.fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_"$TYPE"_"$DOWN"_
	samtools view -@ 5 -Sb -F 260 -q 5 "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam > "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	rm "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam
	# 8.2.2. Assign reads to a feature
	featureCounts -T 1 -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools view -h FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam | awk -F'\t' '$1~/_.{8}$|@SQ|@HD/{print $0}' | samtools view -Sb - > "$ID"_assignedUMI_"$DOWN".bam 
	samtools sort -@ 5 "$ID"_assignedUMI_"$DOWN".bam  > "$ID"_assignedUMI_"$DOWN".sorted.bam
	rm FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam
	samtools index "$ID"_assignedUMI_"$DOWN".sorted.bam
	# 8.2.3. Count UMI
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "$ID"_assignedUMI_"$DOWN".sorted.bam -S FEATURECOUNTS/umi.counts."$DOWN".tsv.gz
	# 8.2.4. Transcript count
	if [ "$DOWN" == "100000" ]
	then
	$salmon quant -i $SALMON_ALL -p 5 -l SF --fldMean 700 --fldSD 200 --fldMax 1100 --minAssignedFrags 1 --dumpEqWeights -r "$TYPE"."$DOWN".R1.fastq --validateMappings -o "$ID"_salmon_UMI_"$DOWN"
	fi
	else
	#### 8.3. INTERNAL ####
	# 8.3.1. Map and count reads: INTERNAL reads
	$STAR --runThreadN 5 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_75" --readFilesIn "$TYPE"."$DOWN".R1.fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_"$TYPE"_"$DOWN"_
	samtools view -@ 5 -Sb -F 260 -q 5 "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam > "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	rm "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam
	# 8.3.2. Assign reads to a feature
	featureCounts -T 1 -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.GENE.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	featureCounts -T 1 -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.EXON.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	# 8.3.3. Transcript count
	if [ "$DOWN" == "100000" ]
	then
	$salmon quant -i $SALMON_ALL -p 5 -l U --fldMean 700 --fldSD 200 --fldMax 1100 --minAssignedFrags 1 --dumpEqWeights -r "$TYPE"."$DOWN".R1.fastq --validateMappings -o "$ID"_salmon_INTERNAL_"$DOWN"
	fi
	fi

	cd "$OUT"/"$ID"/DOWNSAMPLE/

	fi

	done

	cd "$OUT"/"$ID"/

	done

	# 9. ReSQC
	mkdir ReSQC
	# 9.1. GeneBodyCoverage
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_INTERNAL_geneBody.all &
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_UMI_geneBody.all &
	# 9.2. ReSQC Statistics
	# read_distribution.py -i STAR/"$ID"_ALL_Aligned.sortedByCoord.out.bam -r "$ReSQCBED" > ReSQC/"$ID"_ALL_readDistribution.txt &
	read_distribution.py -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_INTERNAL_readDistribution.txt &
	read_distribution.py -i STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_UMI_readDistribution.txt &
	# 9.3. ReSQC various
	read_GC.py -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.filtered.bam-o ReSQC/"$ID"_INTERNAL_GCcontent.txt &
	read_GC.py -i STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_UMI_GCcontent.txt &

	cd "$OUT"
done

# 10. Unload STAR references
STAR --genomeLoad Remove --genomeDir "$STAR_75"
STAR --genomeLoad Remove --genomeDir "$STAR_50"

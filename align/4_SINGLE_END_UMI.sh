#!/bin/bash

# PATH
SAMPLESHEET="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/ID_all.txt"
OUT="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/STAR/"
IN="/home/vincent.hahaut/data_storage/210322_NB551561_0044_AHNTL2AFX2/"

# REFERENCES
#STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/STAR/"
STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/STAR_100bp/" # To use whith 100 bp reads
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/BBDUK/adapters.fa"
GTF_EXONINTRON="/home/vincent.hahaut/data_storage/REFERENCES/hsap_flashseq/REFERENCE/GTF/gencode.v34.exonsANDintrons.gtf"

# BINARIES
STAR="/home/vincent.hahaut/binaries/STAR-2.7.3a/bin/Linux_x86_64/STAR"
BBMAP_filter="/home/vincent.hahaut/data_storage/binaries/bbmap/filterbyname.sh"

mkdir "$OUT"
cd "$OUT"

ulimit -n 100000
STAR --genomeLoad LoadAndExit --genomeDir $STAR_REF

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the variables and FASTQ
	ID="${myArray[0]}"

	# In case multiple folder separated by , are provided.
	FASTQ_R1=$(echo $IN | tr "," "\n" | awk -v ID=$ID '{print $1"/out/"ID"*_R1_*fastq.gz"}')
	#FASTQ_R1="$IN"/out/"$ID".fq.gz

	mkdir "$ID"
	cd "$ID"
	# Pull all FASTQ R1 to the folder:
	cat $FASTQ_R1 > sample.R1.fastq.gz

	echo "$ID"

	# 1. Get UMI
	# Depending on the cell ID - extract the right UMI - Remove adapter/spacer/GGG 
	mkdir FASTQ 

        if echo "$ID" | grep -q "Spacer1"
        then
        echo "TSO-Spacer1"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CAGCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "Spacer2"
        then
        echo "TSO-Spacer2"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>ATAACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "CTAAC"
        then
        echo "TSO-CTAAC"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CTAACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "ATGAC"
        then
        echo "TSO-ATGAC"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>ATGACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "CTGAC"
        then
        echo "TSO-CTGAC"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CTGACGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "AAGCA"
        then
        echo "TSO-AAGCA"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>AAGCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -q "CATCA"
        then
        echo "TSO-CATCA"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*GT)(?P<umi_1>.{8})(?P<discard_2>CATCAGG)(?P<discard_3>G{0,2})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -qE "FStso|Fstso"
	then
	echo "FS-TSO"
   	umi_tools extract --bc-pattern="^(?P<discard_1>.+CAACGCAGAGT){s<=1}(?P<umi_1>.{8})(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	elif echo "$ID" | grep -qE "5_SS3|6_SS3|SS3fwd"
	then
	echo "SS3-TSO"
        umi_tools extract --bc-pattern="^(?P<discard_1>.*ATTGCGCAATG){s<=1}(?P<umi_1>.{8})(?P<discard_2>G{0,4})" --stdin=sample.R1.fastq.gz --stdout=FASTQ/umi.R1.fq --extract-method=regex
	fi

	# UMI-tools has a filtered-out option but does not work here. I think it expect paired ends for this to work.
	# Split UMI - internal using BBMAP
	awk 'NR%4==1{print}' FASTQ/umi.R1.fq | sed 's/_........//g' | sed 's/^@//g' > names.txt
	$BBMAP_filter in=sample.R1.fastq.gz out=FASTQ/internal.R1.fq names=names.txt include=f overwrite=t
	rm names.txt

	cat FASTQ/umi.R1.fq FASTQ/internal.R1.fq > FASTQ/allreads.R1.fq

	# 2. Remove leftover of adapters
  	# t = threads
	# ktrim = right and left
	# rcomp=f only look for the forward sequence not reverse
	# k = min kmer length trim
	# hdist = hamming distance for mismatches
	# mink = min kmer length trim - at edges of fragments
	# hdist2 = for kmer < k, use this hamming distance
	# minlength = minimal read length or discard
	# tbo = trim on where paired reads overlap 
	$BBDUK -Xmx48g in=FASTQ/allreads.R1.fq out=FASTQ/allreads.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 minlength=60 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/allreads.R1.2.fq out=FASTQ/allreads.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 minlength=60 mink=10 hdist=1

  	$BBDUK -Xmx48g in=FASTQ/umi.R1.fq out=FASTQ/umi.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=18 minlength=40 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/umi.R1.2.fq out=FASTQ/umi.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=18 minlength=40 mink=10 hdist=1

  	$BBDUK -Xmx48g in=FASTQ/internal.R1.fq out=FASTQ/internal.R1.2.fq t=32 ktrim=r ref="$BBDUK_REF" k=23 minlength=60 mink=10 hdist=1
	$BBDUK -Xmx48g in=FASTQ/internal.R1.2.fq out=FASTQ/internal.R1.trim.fq t=32 ktrim=l ref="$BBDUK_REF" k=23 minlength=60 mink=10 hdist=1

	rm FASTQ/allreads.R1.2.fq FASTQ/internal.R1.2.fq FASTQ/umi.R1.2.fq

	gzip -c FASTQ/allreads.R1.trim.fq > FASTQ/allreads.R1.fq.gz
	gzip -c FASTQ/umi.R1.trim.fq > FASTQ/umi.R1.fq.gz
	gzip -c FASTQ/internal.R1.trim.fq > FASTQ/internal.R1.fq.gz
	rm FASTQ/allreads.R1.trim.fq FASTQ/umi.R1.trim.fq FASTQ/internal.R1.trim.fq

	# REDEFINED VARIABLES
	MYFASTQ_R1="$OUT"/"$ID"/FASTQ/allreads.R1.fq.gz
	MYFASTQ_UMI_R1="$OUT"/"$ID"/FASTQ/umi.R1.fq.gz
	MYFASTQ_INTERNAL_R1="$OUT"/"$ID"/FASTQ/internal.R1.fq.gz

	# 3. Align the Data
	mkdir STAR

	# To validate trimming results ==> Output no trimmed as well
	# OUTPUT EVERYTHING --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outReadsUnmapped Fastx
	$STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn sample.R1.fastq.gz --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_ALLnoTrim_
	samtools view -@ 30 -Sb -F 260 -q 5 STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.out.bam > STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_ALLnoTrim_Aligned.sortedByCoord.filtered.bam

	for READS in "$MYFASTQ_R1" "$MYFASTQ_UMI_R1" "$MYFASTQ_INTERNAL_R1"
	do

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

	# 3.1. MAP
	$STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn "$READS" --quantMode TranscriptomeSAM --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_"$TYPE"_

	# 3.2. Samtools
	samtools view -@ 30 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam

	samtools view -@ 30 -Sb -F 260 -q 5 STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.out.bam > STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.filtered.bam
	samtools sort -@ 30 STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.filtered.bam -o STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.filtered.sorted.bam
	samtools index STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.filtered.sorted.bam
	rm STAR/"$ID"_"$TYPE"_Aligned.toTranscriptome.filtered.bam

	# 3.3.1. FeatureCounts
	mkdir FEATURECOUNTS
	if [ "$TYPE" != "UMI" ]
	then
	featureCounts -T 15 -t gene -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.GENE.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	featureCounts -T 15 -t exon -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_featureCounts.gencode.EXON.txt STAR/"$ID"_"$TYPE"_Aligned.sortedByCoord.filtered.bam &
	fi

	# 3.3.2. FeatureCounts - BAM file
	# For strand-invasion analysis
	if [ "$TYPE" == "UMI" ]
	then
	# featureCounts -T 30 -t gene -g gene_name -s 0 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_unstranded.Aligned.sorted.bam STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam
	featureCounts -T 30 -t all -g feature -s 0 --largestOverlap --fracOverlap 0.25 -a "$GTF_EXONINTRON" -R BAM -o FEATURECOUNTS/"$ID"_UMI_featureCounts_ExonIntrons STAR/"$ID"_UMI_Aligned.sortedByCoord.filtered.bam 
	cd "$OUT"
	fi

	# 3.4. Downsampling
	mkdir DOWNSAMPLE
	cd DOWNSAMPLE
	for DOWN in 10000 25000 50000 75000 100000 125000 250000
	do
	if [ `zcat "$READS" | grep -c ^@`  -gt "$DOWN" ]
	then
	mkdir DOWN_"$DOWN"
	cd DOWN_"$DOWN"
	seqtk sample -s42 "$READS" "$DOWN" > "$TYPE"."$DOWN".R1.fastq

	# MAP
	$STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn "$TYPE"."$DOWN".R1.fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_"$TYPE"_"$DOWN"_
	samtools view -@ 30 -Sb -F 260 -q 5 "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam > "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools index "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	rm "$ID"_"$TYPE"_"$DOWN"_Aligned.sortedByCoord.out.bam

	# Quantify
	mkdir FEATURECOUNTS
	if [ "$TYPE" = "UMI" ]
	then
	featureCounts -T 30 -t exon -g gene_name -s 1 --fracOverlap 0.25 -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam
	samtools view -h FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam | awk -F'\t' '$1~/_.{8}$|@SQ|@HD/{print $0}' | samtools view -Sb - > "$ID"_assignedUMI_"$DOWN".bam 
	samtools sort -@ 30 "$ID"_assignedUMI_"$DOWN".bam  > "$ID"_assignedUMI_"$DOWN".sorted.bam
	rm FEATURECOUNTS/"$ID"_UMI_"$DOWN"_Aligned.sorted.bam.featureCounts.bam
	samtools index "$ID"_assignedUMI_"$DOWN".sorted.bam
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "$ID"_assignedUMI_"$DOWN".sorted.bam -S FEATURECOUNTS/umi.counts."$DOWN".tsv.gz
	else
	featureCounts -T 15 -t gene -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.GENE.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	featureCounts -T 15 -t exon -g gene_name -a "$GTF" -R BAM -o FEATURECOUNTS/"$ID"_"$TYPE"_"$DOWN".feature.EXON.txt "$ID"_"$TYPE"_"$DOWN"_Aligned.sorted.bam &
	fi

	cd "$OUT"/"$ID"/DOWNSAMPLE/

	fi

	done

	cd "$OUT"/"$ID"/

	done

	# 7. GeneBodyCoverage
	mkdir ReSQC
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.out.bam -o ReSQC/"$ID"_INTERNAL_geneBody.all &
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_ALL_Aligned.sortedByCoord.out.bam -o ReSQC/"$ID"_ALL_geneBody.all &
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_UMI_Aligned.sortedByCoord.out.bam -o ReSQC/"$ID"_UMI_geneBody.all &

	# 8. ReSQC Statistics
	read_distribution.py -i STAR/"$ID"_INTERNAL_Aligned.sortedByCoord.out.bam -r "$ReSQCBED" > ReSQC/"$ID"_INTERNAL_readDistribution.txt &
	read_distribution.py -i STAR/"$ID"_ALL_Aligned.sortedByCoord.out.bam -r "$ReSQCBED" > ReSQC/"$ID"_ALL_readDistribution.txt &
	read_distribution.py -i STAR/"$ID"_UMI_Aligned.sortedByCoord.out.bam -r "$ReSQCBED" > ReSQC/"$ID"_UMI_readDistribution.txt &

	cd "$OUT"
done

STAR --genomeLoad Remove --genomeDir "$STAR_REF"

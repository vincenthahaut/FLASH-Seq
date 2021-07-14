#!/bin/bash

# 1. PATH
# The sample sheet only contains the FASTQ prefix - one per lane
SAMPLESHEET="/home/vincent.hahaut/data_storage/210617_NB551561_0051_AHTYJYAFX2/ID.txt"
OUT="/home/vincent.hahaut/data_storage/210617_NB551561_0051_AHTYJYAFX2/STAR/"
IN="/home/vincent.hahaut/data_storage/210617_NB551561_0051_AHTYJYAFX2/"

# 2. REFERENCES
STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR/"
# STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR_100bp/"	# To use when dealing with 100 bp reads
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
GTF_CODING="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.protein_coding.gtf"
GTF_EXON_CODING="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.exons.gtf"
GTF_INTRON_CODING="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.coding.introns.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BBDUK/adapters.fa"

# 3. BINARIES
STAR="/home/vincent.hahaut/binaries/STAR-2.7.3a/bin/Linux_x86_64/STAR"
BBDUK="/home/vincent.hahaut/data_storage/binaries/bbmap/bbduk.sh"
BBSPLIT="/home/vincent.hahaut/data_storage/binaries/bbmap/bbsplit.sh"
BBMAP_filter="/home/vincent.hahaut/data_storage/binaries/bbmap/filterbyname.sh"

# 4. Run the pipeline for every sample
cd "$IN"
mkdir "$OUT"
cd "$OUT"

ulimit -n 100000
STAR --genomeLoad LoadAndExit --genomeDir $STAR_REF --runThreadN 24

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the CELL ID and FASTQ
	ID="${myArray[0]}"
	FASTQ=$(echo $IN/out/$ID*R1*fastq.gz | sed 's/ /,/g')
	echo "$ID"
	mkdir "$ID"
	cd "$ID"

	# 1. Save the FASTQ
	cp $(echo "$FASTQ" | sed 's/,/ /g') .
	cat *fastq.gz > sample.fastq.gz

	# 2. Trim Reads
	if echo "$ID" | grep -q "SS3"
	then
	#### Remove UMI from 5'UMI reads from SS3 data and add them to the read name
	# t = threads
	# ktrim = right and left
	# rcomp=f only look for the forward sequence not reverse
	# k = min kmer length trim
	# hdist = hamming distance for mismatches
	# mink = min kmer length trim - at edges of fragments
	# hdist2 = for kmer < k, use this hamming distance
	# minlength = minimal read length or discard
	# tbo = trim on where paired reads overlap 
	# 2.1. Remove UMI
	# Could leave some reads with .*ATTGCGCAATG{umi} but these are rare and also left by zUMIs
    	umi_tools extract --bc-pattern="^(?P<discard_1>ATTGCGCAATG){s<=2}(?P<umi_1>.{8})" --stdin=sample.fastq.gz --stdout=umi.fastq --extract-method=regex
	# 2.2. Subset non-UMI reads and append them to trimmed umi.fastq
	awk 'NR%4==1{print}' umi.fastq | sed 's/_........//g' | sed 's/^@//g' > names.txt
	$BBMAP_filter in=sample.fastq.gz out=non_umi.fastq names=names.txt include=f overwrite=t
  	cat non_umi.fastq umi.fastq > preprocessed.fq
    	# 2.3. Trim Reads
    	# Remove F/R sequencing adapters
    	$BBDUK -Xmx48g in=preprocessed.fq out=cleaned.polyA.fastq t=32 ktrim=r rcomp=f literal=AAAAAAA hdist=0 k=6
    	$BBDUK -Xmx48g in=cleaned.polyA.fastq out=cleaned.left.fastq t=32 ktrim=l ref="$BBDUK_REF" k=23 minlength=50 mink=7 hdist=1 hdist2=0 tbo
    	$BBDUK -Xmx48g in=cleaned.left.fastq out=cleaned.fastq t=32 ktrim=r ref="$BBDUK_REF" k=23 minlength=50 mink=7 hdist=1 hdist2=0 tbo
    	else
    	#### Takara, SS2, FS
    	# 2.3. Trim Reads
    	$BBDUK -Xmx48g in=sample.fastq.gz out=cleaned.left.fastq t=32 ktrim=l ref="$BBDUK_REF" k=23 minlength=60 mink=7 hdist=1 hdist2=0 tbo
    	$BBDUK -Xmx48g in=cleaned.left.fastq out=cleaned.fastq t=32 ktrim=r ref="$BBDUK_REF" k=23 minlength=60 mink=7 hdist=1 hdist2=0 tbo
	fi
	gzip -c cleaned.fastq > cleaned.fastq.gz

	mkdir FASTQ
	MYFASTQ="$OUT"/"$ID"/cleaned.fastq.gz
	scp "$MYFASTQ" FASTQ/

	# 3. Align the Data
	# TO OUTPUT EVERYTHING: --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outReadsUnmapped Fastx
	mkdir STAR
	$STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn "$MYFASTQ" --readFilesCommand zcat --quantMode TranscriptomeSAM --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_

	# 4. Samtools
	samtools view -@ 30 -Sb -F 260 STAR/"$ID"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_Aligned.sortedByCoord.filtered.bam

	samtools view -@ 30 -Sb -F 260 STAR/"$ID"_Aligned.toTranscriptome.out.bam > STAR/"$ID"_Aligned.toTranscriptome.filtered.bam
	samtools sort -@ 30 STAR/"$ID"_Aligned.toTranscriptome.filtered.bam -o STAR/"$ID"_Aligned.toTranscriptome.filtered.sorted.bam
	samtools index STAR/"$ID"_Aligned.toTranscriptome.filtered.sorted.bam
	rm STAR/"$ID"_Aligned.toTranscriptome.filtered.bam

	# 5. Read Count using FeatureCounts
	# -T = Threads
	# -t = Summarise exon results mapping ..
	# -g = .. to gene IDs
	# -a = reference
	# -o = OUTPUT
	mkdir FEATURECOUNTS
	# On EXON - ALL GENES
	featureCounts -T 1 -t exon -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &
	# On GENE - ALL GENES
	featureCounts -T 1 -t gene -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.GENE.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &
	# On EXON - CODING
	featureCounts -T 1 -t exon -g gene_name -a "$GTF_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &
	# On GENE - CODING
	featureCounts -T 1 -t gene -g gene_name -a "$GTF_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.GENE.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &

	# 6. Downsampling
	mkdir DOWNSAMPLE
	cd DOWNSAMPLE
	for DOWN in 50000 100000 250000 500000 1000000
	do
	if [ `zcat "$MYFASTQ" | grep -c ^@`  -gt "$DOWN" ]
	then
	# If multiple downsamplings were needed
	for TIMES in 1
	do
	mkdir DOWN_"$DOWN"_"$TIMES"
	cd DOWN_"$DOWN"_"$TIMES"
	if [ $TIMES -eq 1 ]
	then
	seqtk sample -s42 "$MYFASTQ" "$DOWN" > downsampled."$TIMES".fastq
	fi
	#
	# 6.2. STAR Mapping
	mkdir FEATURECOUNTS
	$STAR --runThreadN 30 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --quantMode TranscriptomeSAM --readFilesIn downsampled."$TIMES".fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_downsampled_"$DOWN"_"$TIMES"_
	samtools view -@ 30 -Sb -F 260 "$ID"_downsampled_"$DOWN"_"$TIMES"_Aligned.sortedByCoord.out.bam > "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam
	samtools index "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam

	# 6.3. FeatureCounts Quantification
	# Many of the following lines are not directly used in the paper but were used to explore the results.
	# On EXON - ALL GENES
	featureCounts -T 1 -t exon -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On GENE - ALL GENES
	featureCounts -T 1 -t gene -g gene_name -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.GENE."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On EXON - CODING
	featureCounts -T 1 -t exon -g gene_name -a "$GTF_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On GENE - CODING
	featureCounts -T 1 -t gene -g gene_name -a "$GTF_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.GENE."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On EXON - CODING
	featureCounts -T 1 -t exon -g gene_id -a "$GTF_EXON_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.exonOnly."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On INTRON - CODING
	# Read need to be fully enclosed in a intron to be called intronic
	featureCounts -T 1 -t intron -g gene_id -a "$GTF_INTRON_CODING" --fracOverlap 1 -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.intronOnly."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On EXON - CODING
	featureCounts -T 1 -t exon -g exon_id -a "$GTF_EXON_CODING" -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.exonOnly.SummarizedToExonID."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On INTRON - CODING
	featureCounts -T 1 -t intron -g intron_id -a "$GTF_INTRON_CODING" --fracOverlap 1 -o FEATURECOUNTS/"$ID"_ReadCount_coding.featureCounts.gencode.intronOnly.SummarizedToIntronID."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &

	# 6.4. Clean-up
	rm *_Aligned.filtered.bam *_Aligned.sortedByCoord.out.bam *_Log.out *Log.progress.out *_SJ.out.tabs
	cd ..
	done
	fi
	done
	cd ..

	# 7. GeneBodyCoverage
	mkdir ReSQC
	geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_geneBody.all &

	# 8. ReSQC Statistics
	read_distribution.py -i STAR/"$ID"_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_readDistribution.txt &

	# 9. Clean-up
	rm *fastq*
	rm STAR/"$ID"_Aligned.sortedByCoord.filtered.bam.tmp
	rm FASTQC/*zip
	rm log.txt

	cd "$OUT"
done

STAR --genomeLoad Remove --genomeDir "$STAR_REF"

rm ./*/DOWNSAMPLE/*/down*fastq

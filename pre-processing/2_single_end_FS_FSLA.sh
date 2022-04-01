#!/bin/bash

ulimit -n 100000

# 0. PATH
# The sample sheet only contains the FASTQ prefix - one per lane
SAMPLESHEET="/home/vincent.hahaut/data_storage/210122_NB551561_0039_AHH5VCBGXH/xad"
OUT="/home/vincent.hahaut/data_storage/210122_NB551561_0039_AHH5VCBGXH/STAR/"
IN="/home/vincent.hahaut/data_storage/210122_NB551561_0039_AHH5VCBGXH/"

# 1. REFERENCES
STAR_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/STAR/"
GTF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.annotation.gtf"
GTF_EXONINTRON="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/GTF/gencode.v34.primary_assembly.exonsANDintrons.gtf"
ReSQCBED="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BED/hg38_UCSC_gencodeV34.ReSQCBED.bed"
BBDUK_REF="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/REFERENCE/BBDUK/adapters.fa"
SALMON_ALL="/home/vincent.hahaut/data_storage/REFERENCES/hsap_fastsmart/hsap_fastsmart/salmon_precomputed/salmon_sa_index/"

# 3. BINARIES
STAR="/home/vincent.hahaut/data_storage/binaries/STAR/bin/Linux_x86_64/STAR"
BBDUK="/home/vincent.hahaut/data_storage/binaries/bbmap/bbduk.sh"
BBMAP_filter="/home/vincent.hahaut/data_storage/binaries/bbmap/filterbyname.sh"
TRIM="/home/vincent.hahaut/data_storage/binaries/Trimmomatic-0.39/trimmomatic-0.39.jar"
salmon=/home/vincent.hahaut/data_storage/binaries/salmon-1.5.2_linux_x86_64/bin/salmon

# 4. Run the pipeline for every sample
cd "$IN"
mkdir "$OUT"
cd "$OUT"

STAR --genomeLoad LoadAndExit --genomeDir $STAR_REF --runThreadN 24

cat "$SAMPLESHEET" |
while IFS=$'\t' read -r -a myArray; do

	# 0. Get the CELL ID
	ID="${myArray[0]}"
	echo "$ID"
	mkdir "$ID"
	cd "$ID"

	# 1. Get the FASTQ
	FASTQ=`ls $IN/out/${ID}*R1*fastq.gz`
	cp $FASTQ .
	cat *fastq.gz > sample.fastq.gz

	# 2. Trim Reads
	if echo "$ID" | grep -q "SS3"
	then
		#### Remove UMI from 5'UMI reads from SS3 data and add them to the read name
		# This pipeline assumes that the UMI will not be used later. Refers to SINGLE_END_UMI.sh for more details on UMI
		# 2.1. Remove UMI
    	umi_tools extract --bc-pattern="^(?P<discard_1>ATTGCGCAATG){s<=1}(?P<umi_1>.{8})" --stdin=sample.fastq.gz --stdout=umi.fastq --extract-method=regex
		# 2.2. Subset non-UMI reads and append them to trimmed umi.fastq
		awk 'NR%4==1{print}' umi.fastq | sed 's/_........//g' | sed 's/^@//g' > names.txt
		$BBMAP_filter in=sample.fastq.gz out=non_umi.fastq names=names.txt include=f overwrite=t
  		cat non_umi.fastq umi.fastq > preprocessed.fq
    	# 2.3. Trim Reads
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
    	$BBDUK -Xmx48g in=preprocessed.fq out=cleaned.left.fastq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    	$BBDUK -Xmx48g in=cleaned.left.fastq out=cleaned.fastq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    	else
    	#### Takara, SS2, FS
    	# 2.3. Trim Reads
    	$BBDUK -Xmx48g in=sample.fastq.gz out=cleaned.left.fastq t=32 ktrim=l ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
    	$BBDUK -Xmx48g in=cleaned.left.fastq out=cleaned.fastq t=32 ktrim=r ref="$BBDUK_REF" k=23 mink=7 hdist=1 hdist2=0 tbo
		fi

	# 3. Trim Reads to 75bp if longer
	# Just used for a fair comparison
	java -jar $TRIM SE -threads 10 -phred33 cleaned.fastq cleaned2.fastq CROP:75 MINLEN:45
	mv cleaned2.fastq cleaned.fastq

	# 4. Save the Cleaned FASTQ file
	mkdir FASTQ 
	gzip -c cleaned.fastq > cleaned.fastq.gz
	MYFASTQ="$OUT"/"$ID"/cleaned.fastq.gz
	scp "$MYFASTQ" FASTQ/

	# 5. Align the Data
	# TO OUTPUT EVERYTHING: --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outReadsUnmapped Fastx
	mkdir STAR
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn "$MYFASTQ" --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR/"$ID"_

	# 6. Samtools
	# Filter out secundary alignments and unmapped reads
	samtools view -@ 10 -Sb -F 260 STAR/"$ID"_Aligned.sortedByCoord.out.bam > STAR/"$ID"_Aligned.sortedByCoord.filtered.bam
	samtools index STAR/"$ID"_Aligned.sortedByCoord.filtered.bam

	# 6. Read Count using FeatureCounts
	# -T = Threads
	# -t = Summarise exon results mapping ..
	# -g = .. to gene IDs
	# -a = reference
	# -o = OUTPUT
	mkdir FEATURECOUNTS
	# On EXON - ALL GENES
	featureCounts -T 1 -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &
	# On GENE - ALL GENES
	featureCounts -T 1 -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.GENE.txt STAR/"$ID"_Aligned.sortedByCoord.filtered.bam &

	# 7. Downsampling
	mkdir DOWNSAMPLE
	cd DOWNSAMPLE
#	for DOWN in 5000 10000 20000 40000 50000 75000 100000 125000 250000 375000 500000
	for DOWN in 100000 250000
	do
	if [ `zcat "$MYFASTQ" | grep -c ^@`  -gt "$DOWN" ]
	then
	# 7.1. Downsample FASTQ
	# If multiple iterations of the downsampling were needed (not the case)
	for TIMES in 1
	do
	mkdir DOWN_"$DOWN"_"$TIMES"
	cd DOWN_"$DOWN"_"$TIMES"
	if [ $TIMES -eq 1 ]
	then
	seqtk sample -s42 "$MYFASTQ" "$DOWN" > downsampled."$TIMES".fastq
	fi
	# 
	# 7.2. STAR Mapping
	mkdir FEATURECOUNTS
	$STAR --runThreadN 10 --limitBAMsortRAM 20000000000 --genomeLoad LoadAndKeep --genomeDir "$STAR_REF" --readFilesIn downsampled."$TIMES".fastq --readFilesCommand cat --limitSjdbInsertNsj 2000000 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$ID"_downsampled_"$DOWN"_"$TIMES"_
	samtools view -@ 10 -Sb -F 260 "$ID"_downsampled_"$DOWN"_"$TIMES"_Aligned.sortedByCoord.out.bam > "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam
	samtools index "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam
	
	# 7.3. FeatureCounts Quantification
	# On EXON - ALL GENES
	featureCounts -T 1 -t exon -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
	# On GENE - ALL GENES
	featureCounts -T 1 -t gene -g gene_name --fracOverlap 0.25 -a "$GTF" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.GENE."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &

	# 7.4. Transcript quantification using Salmon
	if [ $DOWN -eq 250000 ]
	then
	# 6.2. Map and Quantify using salmon
	# -i = index
	# -p = threads
	# -l U (libtypes) = unstranded single-ends reads expected
	# -r = single-end input (vs -1 and -2 for paired)
	# --validateMappings = improve both the sensitivity and specificity of mapping
	# -o = output prefix
	# --fldMean = What is the fragment length mean ? (From bioanalyzer)
	# --fldSD = see previous
	# --fldMax = Max size fragments
	# --minAssignedFrags = minimal number match to quantify
	# Not used --geneMap = if GTF specified, summary transcripts at the gene level  
 	$salmon quant -i $SALMON_ALL -p 4 -l U --fldMean 500 --fldSD 300 --fldMax 1100 --minAssignedFrags 1 --dumpEqWeights -r downsampled."$TIMES".fastq --validateMappings -o salmon_ALL_TRANSCRIPTS_"$TIMES"
	
	# 7.5. Intron - Exon quantification to explore the results
	featureCounts -T 1 -t all -g feature_name -s 0 --fracOverlap 0.25 --largestOverlap -a "$GTF_EXONINTRON" -o FEATURECOUNTS/"$ID"_ReadCount.featureCounts.gencode.EXONINTRON."$DOWN"."$TIMES".txt "$ID"_downsampled_"$DOWN"_Aligned.sortedByCoord.filtered."$TIMES".bam &
    fi
    
	# 7.6. Clean-up
	rm *_Log.out *Log.progress.out
	cd ..
	done
	fi
	done
	cd ..

	# 8. ReSQC
	mkdir ReSQC

	# 8.1. GeneBodyCoverage
	# Slow down the script ++ - turn off if not needed
	# geneBody_coverage.py -r "$ReSQCBED" -i STAR/"$ID"_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_geneBody.all &

	# 8.2. ReSQC Statistics
	read_distribution.py -i STAR/"$ID"_Aligned.sortedByCoord.filtered.bam -r "$ReSQCBED" > ReSQC/"$ID"_readDistribution.txt &

	# 8.3. ReSQC various
	read_GC.py -i STAR/"$ID"_Aligned.sortedByCoord.filtered.bam -o ReSQC/"$ID"_GCcontent.txt &

	# 10. Clean-up
	rm *fastq*
	rm STAR/"$ID"_Aligned.sortedByCoord.filtered.bam.tmp
	rm FASTQC/*zip
	rm log.txt

	cd "$OUT"
done

# 11. Unload STAR Reference
STAR --genomeLoad Remove --genomeDir "$STAR_REF"

# 12. Final Clean-up
rm ./*/DOWNSAMPLE/*/down*fastq

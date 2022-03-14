#!/bin/bash


# 0. PATH
# Only contains R1 Fastq files from selected cells
FASTQ_PATH="/home/vincent.hahaut/data_storage/TRUST4/FASTQ/"

# 1. BINARIES
TRIM="/home/vincent.hahaut/data_storage/binaries/Trimmomatic-0.39/trimmomatic-0.39.jar"

# 2. DOWNSAMPLING
cd "$FASTQ_PATH"
for READS in *
do
echo $READS

# 2.1. Trim reads
# Some of them were sequenced to 100bp R1
java -jar $TRIM SE -threads 10 -phred33 "$READS" trimmed.fq CROP:75 MINLEN:70

# 2.2. Downsample
for DOWN in 5000 10000 20000 40000 50000 75000 100000 125000 250000 375000 500000
do
mkdir DOWN_"$DOWN"
if [ `zcat "$READS" | grep -c ^@`  -gt "$DOWN" ]
then
seqtk sample -s42 trimmed.fq "$DOWN" > DOWN_"$DOWN"/"$READS"."$DOWN".fastq
fi
done
rm trimmed.fq

done


# 3. Run TRUST4 on each downsampled set
cd /home/vincent.hahaut/data_storage/TRUST4/

# 3.1. Get the list of downsampled fastq
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_5000" > trust4.5000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_10000" > trust4.10000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_20000" > trust4.20000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_40000" > trust4.40000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_50000" > trust4.50000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_75000" > trust4.75000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_100000" > trust4.100000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_125000" > trust4.125000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_250000" > trust4.250000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_375000" > trust4.375000.txt
ls "/home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_500000" > trust4.500000.txt

# 3.2. Run TRUST4-smart-seq for each downsampling
cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_5000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.5000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_5K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_10000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.10000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_10K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_20000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.20000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_20K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_40000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.40000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_40K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_50000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.50000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_50K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_75000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.75000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_75K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_100000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.100000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_100K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_125000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.125000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_125K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_250000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.250000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_250K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_375000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.375000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_375K

cd /home/vincent.hahaut/data_storage/TRUST4/FASTQ/DOWN_500000
perl /home/vincent.hahaut/data_storage/binaries/TRUST4/trust-smartseq.pl -1 /home/vincent.hahaut/data_storage/TRUST4/trust4.500000.txt -f /home/vincent.hahaut/data_storage/binaries/TRUST4/hg38_bcrtcr.fa --ref /home/vincent.hahaut/data_storage/binaries/TRUST4/human_IMGT+C.fa -t 8 -o /home/vincent.hahaut/data_storage/TRUST4/TRUST_500K

#!/bin/bash
#Aligns Wamaitha_2015 Gata6 ChIPseq data to mm10

path=$1
fastq_dir=$2
work_dir=$3
ebwt=$4

cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd ${fastq_dir}

for f in $(ls *.fastq | rev | cut -c 9- | rev | uniq)
do
	cd ${data_dir}

	echo -e "Trimming $f with trim_galore"
	prun python3 trim_galore --paired --illumina ${fastq_dir}$f\_1.fastq ${fastq_dir}$f\_2.fastq >> $f\.trimmingStats.txt 2>&1

	echo -e "Mapping $f with bowtie2"
	bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 2000 -x $ebwt -1 $f\_1_val_1.fq -2 $f\_2_val_2.fq -S $f\.sam >> $f\.mappingStats.txt 2>&1
	rm $f\_1_val_1.fq
	rm $f\_2_val_2.fq

	echo -e "Creating BAM files for $f and filtering for paired and mapped reads"
	samtools view -b -h -f 2 -q 20 $f\.sam > $f\.bam
	rm $f\.sam

	echo -e "Sorting BAM files for $f"
	samtools sort -m 1G $f\.bam -T $f\_sorted -o $f\_sorted.bam

	echo -e "Removing Blacklisted regions from $f"
	bedtools intersect -v -a $f\_sorted.bam -b  ${path}files/mm10.bl.bed > $f\_sorted_blacklisted.bam

	echo -e "Removing duplicates from $f using PICARD\n"
	java -jar $PICARD MarkDuplicates INPUT=$f\_sorted_blacklisted.bam \
			OUTPUT=${bam_dir}$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT \
			REMOVE_DUPLICATES=TRUE
	samtools index ${bam_dir}$f\_dedup.bam
done

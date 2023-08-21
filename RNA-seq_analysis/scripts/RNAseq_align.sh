#Pipeline for analysis of RNA-seq data (TX timecourse 2016)
path=$1
fastq_dir=$2
work_dir=$3
genome_dir=$4

cd $work_dir

mkdir -p data
data_dir=${work_dir}data'/'

mkdir -p final_bam
bam_dir=${work_dir}final_bam'/'

cd $fastq_dir

for f in $(ls *.fastq | rev | cut -c 14- | rev | uniq)
do
	cd ${data_dir}

  echo -e "Mapping $f with STAR"
  STAR --genomeDir $genome_dir --readFilesIn ${fastq_dir}$f\_R1_001.fastq ${fastq_dir}$f\_R2_001.fastq \
				--outSAMtype BAM Unsorted --outFileNamePrefix $f\_ --outSAMattributes NH HI NM MD

  echo -e "Filtering $f for properly paired reads"
  samtools view -q 7 -f 3 -b $f\_Aligned.out.bam > $f\_filtered.bam

  echo -e "Sorting BAM files for $f"
  samtools sort -m 1G $f\_filtered.bam -T $f\_sorted -o ${bam_dir}$f\_sorted.bam
  samtools index ${bam_dir}$f\_sorted.bam

done

#!/bin/bash
# Merges files ending on _dedup.bam

work_dir=$1
bam_dir=$2 # Contains BAM files
cd $work_dir

mkdir -p merged_bam
merged_bam_dir=${work_dir}merged_bam'/'

cd $bam_dir

for f in $(ls *_r[1-2]_dedup.bam | rev | cut -c 14- | rev | uniq)
do
  echo -e "Merging BAM files for $f"
  samtools merge ${merged_bam_dir}$f\_merged.bam $f\_r1_dedup.bam \
      $f\_r2_dedup.bam
  samtools index ${merged_bam_dir}$f\_merged.bam
done

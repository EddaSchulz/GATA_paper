#!/bin/bash
# Creates BIGWIG tracks for unspliced NGS data (ATAC-seq, CUT&Tag, ChIP-seq, TT-seq)


work_dir=$1
bam_dir=$2 # Contains BAM files

cd $work_dir

mkdir -p bigwig
bigwig_dir=${work_dir}bigwig'/'

cd $bam_dir

for f in $(ls *.bam | rev | cut -c 5- | rev | uniq)
do
  echo -e "Creating BIGWIG tracks for $f"
  prun python3 bamCoverage -b ${bam_dir}$f\.bam -o ${bigwig_dir}$f\.bw -bs 10 -e \
    --normalizeUsing CPM -ignore chrX chrY
done

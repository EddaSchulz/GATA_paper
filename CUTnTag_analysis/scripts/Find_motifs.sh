#!/bin/bash
path=$1
peaks=$2
wd=$3


mkdir -p ${wd}fimo_data'/'

data_dir=${wd}fimo_data/
output_dir=/project/agsgpa/Till/CutnTag/Liat_Gata/output_final/gata_peaks_v2/output/
motif_dir=/project/agsgpa/Till/CutnTag/Liat_Gata/output_final/meme/jaspar_motifs/
peak_dir=/project/agsgpa/Till/CutnTag/Liat_Gata/output_final/merged_peaks/

cd $peaks

for f in $(ls *Gata* | rev | cut -c 5- | rev | uniq)
do
  g=$(echo -e "$f" | grep -o "Gata[0-9]")
  bedtools getfasta -fi ${path}files/mm10.fa -bed $f\.bed > ${data_dir}$f\.fa
  echo -e "fimo for $f"
  fimo --oc ${data_dir}$f\_fimo --max-stored-scores 5000000 --thresh 0.001 ${path}files/$g\_motif.meme ${data_dir}$f\.fa
  awk -v OFS='\t' '{print $3}' ${data_dir}$f\_fimo/fimo.tsv | grep "chr" | \
      sed 's/[:-]/\t/g' > ${data_dir}$f\_gata_motifs_unsorted.bed
  bedtools sort -i ${data_dir}$f\_gata_motifs_unsorted.bed > ${data_dir}$f\_gata_motifs.bed

  echo -e "Splits peak files depending on GATA_motifs for $f"
  bedtools intersect -u -a $f\.bed -b ${data_dir}$f\_gata_motifs.bed \
      | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4=".", $5="1000", $6=".", $7=$2, $8=$3, $9="234,010,142"}' \
      > ${data_dir}$f\_gata_motif_pos.bed
  bedtools intersect -v -a $f\.bed -b ${data_dir}$f\_gata_motifs.bed \
      | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4=".", $5="1000", $6=".", $7=$2, $8=$3, $9="150,150,150"}' \
      > ${data_dir}$f\_gata_motif_neg.bed

  echo -e "Combines differential and consensus peaks for $f"
  cat ${data_dir}$f\_gata_motif_pos.bed ${data_dir}$f\_gata_motif_neg.bed > ${data_dir}$f\_rgb_motif_unsorted.bed
  bedtools sort -i ${data_dir}$f\_rgb_motif_unsorted.bed > ${wd}$f\_rgb_motif.bed
done

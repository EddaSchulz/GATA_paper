#!/bin/bash

bam_dir=$1
work_dir=$2

cd $bam_dir

bam_files=$(ls *.bam | tr "\n" " ")

prun python3 multiBamSummary bins --smartLabels -e -b $bam_files -o ${work_dir}Gata_multiBam.npz
prun python3 plotCorrelation -in ${work_dir}Gata_multiBam.npz -c pearson -p heatmap -o ${work_dir}Gata_plotCorrelation.pdf \
  --outFileCorMatrix ${work_dir}Gata_CnT_pearson.txt

tail -n +2 Gata_CnT_pearson.txt > Gata_CnT_pearson.txt.tmp && mv Gata_CnT_pearson.txt.tmp Gata_CnT_pearson.txt
sed "s/'//g" Gata_CnT_pearson.txt > Gata_CnT_pearson.txt.tmp && mv Gata_CnT_pearson.txt.tmp Gata_CnT_pearson.txt

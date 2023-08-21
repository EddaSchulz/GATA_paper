#!/bin/bash
path=$1
wd=$2


mkdir -p ${wd}ame_data'/'
data_dir=${wd}ame_data/

cd $wd

for f in $(ls *_ame_input.bed | rev | cut -c 15- | rev | uniq)
do
  echo -e "getfasta for $f"
  bedtools getfasta -fi ${path}files/mm10.fa -bed ${wd}$f\_ame_input.bed > ${data_dir}$f\.fa
  echo -e "ame for $f"
  ame -oc ${data_dir}$f --shuffle-- ${data_dir}$f\.fa ${path}files/JASAPAR2020_vertebrate_non_redundant_PMF_MEME.txt
done

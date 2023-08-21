#!/bin/bash
# Script generates a count table from FASTQ files of the CRISPRa screen

fastq_dir=$1 # FASTQ files are found at GSE194018
library_file=$2 # Path to sgRNA library file (Can be found at /CRISPRa_screen_analysis/files/library_file.txt)
control_file=$3 # Path to list of non-targeting controls (Can be found at /CRISPR_screen_analysis/files/nt_file.txt)
output_dir=$4

cd $fastq_dir

fastq_files=$(find -type f -regex ".*.fastq" -printf "%f ")
names=$(echo -e "$fastq_files" | sed "s/.fastq//g" | sed "s/ /,/g" |  sed "s/.$//g")

mageck count -l $library_file --fastq $fastq_files --norm-method control --control-sgrna $control_file  --pdf-report --sample-label $names \
    --output-prefix ${output_dir}CRISPRa_screen

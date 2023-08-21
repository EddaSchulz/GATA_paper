#!/bin/bash
# Script analyzes count tables generated from MAGeCK count using MAGeCK mle

count_table=$1 # Path to duplicated count table produced with MAGECK_count.sh
design_matrix=$2 # Path to design matrix at /CRISPRa_screen_analysis/files/mle_matrix.txt
control_file=$3 # Path to list of non-targeting at /CRISPR_screen_analysis/files/NTCs_index.txt
output_dir=$4

mageck mle -k $count_table -d $design_matrix --control-sgrna $control_file --norm-method control --output-prefix ${output_dir}CRISPRa_screen

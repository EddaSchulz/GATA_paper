#!/bin/bash
# Master script to run MAGeCK count and mle on a transcript level
# Selects top isoform to re-run MAGeCK on a gene level
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Uses MAGeCK to align CRISPR screen data and perform read counting and statistical analysis."
   echo
   echo "Syntax: ./master1_MAGeCK_analysis.sh [-f|d|p|h]"
   echo "options:"
   echo "f     Provide directory containing FASTQ files (CRISPR Screen). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gata_paper/CRISPRa_screen_analysis/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -f)
            fastq_dir="$2"
            ;;
        -d)
            work_dir="$2"
            ;;
        -p)
            path="$2"
            ;;
        -h)
            help
            exit 0
            ;;
        *)
            echo "Unknown or badly placed parameter '$1'." 1>&2
            exit 1
            ;;
    esac
}

while [[ "$#" -ge 1 ]]; do
    parse_args "$1" "$2"
    shift; shift
done

if [[ $path == '' ]]
then
	echo -e "Please provide the path to /Gata_paper/CRISPRa_screen_analysis/ with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing CRISPRa screen FASTA files with -f"
  exit 1
fi


fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

mkdir -p ${work_dir}mageck_output_raw
mkdir -p ${work_dir}mageck_output_dup
mkdir -p ${work_dir}mageck_output_final


# Aligns reads and performs counting using MAGeCK
echo -e "Counts reads using MAGeCK count"
${path}scripts/MAGeCK_count.sh $fastq_dir ${path}files/Library_file_dup.txt ${path}files/NTCs_index.txt ${work_dir}mageck_output_raw'/'

raw_counts=${work_dir}mageck_output_raw/CRISPRa_screen.count.txt

# Duplicates sgRNAs targeting multiple isoforms
echo -e "Duplicating multi-targeting sgRNAs in count table"
Rscript ${path}scripts/Dup_count_table.R ${work_dir} ${path}

dup_counts=${work_dir}mageck_output_dup/CRISPRa_screen.dup_count.txt

# Performs statistical analysis on duplicated counts using mageck mle
echo -e "Running MAGeCK mle for duplicated counts"
${path}scripts/MAGeCK_mle.sh $dup_counts ${path}files/mle_matrix.txt ${path}files/NTCs_index.txt ${work_dir}mageck_output_dup'/'

# Reducing count tables by top isoform, to re-run mle on a gene level
echo -e "Reducing count tables to single isoforms/gene"
Rscript ${path}scripts/Red_count_table.R ${work_dir} ${path}

final_counts=${work_dir}mageck_output_final/CRISPRa_screen.final_count.txt

# Performs statistical analysis on duplicated counts using mageck mle on the gene level
echo -e "Running MAGeCK mle for final counts"
${path}scripts/MAGeCK_mle.sh $final_counts ${path}files/mle_matrix.txt ${path}files/NTCs_index.txt ${work_dir}mageck_output_final'/'

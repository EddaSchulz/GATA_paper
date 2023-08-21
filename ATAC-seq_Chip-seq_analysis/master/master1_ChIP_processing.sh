#!/bin/bash
# Master script to align Wamaitha 2015 GATA6 Chip-seq
# Calls peaks and creates bigwig tracks
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds mm10 genome and aligns Wamaitha 2015 GATA6 Chip-seq data using Bowtie2. Calls peaks and creates browser tracks."
   echo
   echo "Syntax: ./master1_ChIP_processing.sh [-f|d|p|h]"
   echo "options:"
   echo "f     Provide directory containing FASTQ files (Wamaitha 2015 GATA6 Chip-seq). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gata_paper/ATAC-seq_Chip-seq_analysis/. [mandatory]"
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
	echo -e "Please provide the path to /Gata_paper/ATAC-seq_Chip-seq_analysis/ with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing Wamaitha 2015 GATA6 Chip-seq FASTQ files with -f"
  exit 1
fi


fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/mm10.fa

# Build genome index for bowtie2 from N_masked mm10 genome
echo -e "Build mm10 genome for Bowtie2"
${path}scripts/build_bowtie2.sh $path $work_dir $genome

ebwt=${work_dir}genome/mm10

# Aligns FASTQ files and performs filtering steps
echo -e "Align ChIPseq data"
${path}scripts/ChIP_align.sh $path $fastq_dir $work_dir $ebwt

bam_dir=${work_dir}final_bam'/'

# Calls peaks for replicates and merges them
echo -e "Call peaks from ChIP"
${path}scripts/ChIP_peaks.sh $work_dir $bam_dir

# Merges replicate BAM files
echo -e "Merge replicate BAM files"
${path}scripts/merge_bam.sh $work_dir $bam_dir

merged_bam_dir=${work_dir}merged_bam'/'

# Generates coverage tracks
echo -e "Generate coverage tracks"
${path}scripts/create_bigwig.sh $work_dir $merged_bam_dir

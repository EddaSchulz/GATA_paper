#!/bin/bash
# Master script to align TX timecourse RNA-seq data
# Calculates gene expression (TPM)
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Builds mm10 genome and aligns TX timecourse RNA-seq data using STAR. Calculates TPM using featureCounts."
   echo
   echo "Syntax: ./master1_ChIP_processing.sh [-f|d|p|h]"
   echo "options:"
   echo "f     Provide directory containing FASTQ files (GSE151009, GSE194018). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gata_paper/RNA-seq_analysis/. [mandatory]"
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
	echo -e "Please provide the path to /Gata_paper/RNA-seq_analysis/ with -p"
  exit 1
fi

if [[ $fastq_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing TX timecourse RNA-seq FASTQ files (GSE151009, GSE194018) with -f"
  exit 1
fi


fastq_dir=$(realpath $fastq_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

genome=${path}files/N_masked_B6_Cast.fa
gtf=${path}files/GENCODE_vM25_plus_Xert.gtf

# Build genome index for STAR from N_masked mm10 genome
echo -e "Build mm10 genome for STAR"
${path}scripts/STAR_genomeGenerate.sh $path $work_dir $genome $gtg

ebwt=${work_dir}genome/

# Aligns FASTQ files with STAR
echo -e "Align RNA-seq data"
${path}scripts/RNAseq_align.sh $path $fastq_dir $work_dir $ebwt

bam_dir=${work_dir}final_bam'/'

# Calculate TPM
echo -e "Generate coverage tracks"
Rscript ${path}scripts/TPM_calc.R $bam_dir $work_dir

#!/bin/bash
# Master script to perform qc and analysis for GATA CUT&Tag data
input_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Produces plots for CUT&Tag analysis including correlation and motif enrichment."
   echo
   echo "Syntax: ./master2_CUTnTag_analysis.sh [-b|d|p|h]"
   echo "options:"
   echo "b     Provide directory containing output from master1_CUTnTag_processing.sh (GATA CUT&Tag). [mandatory]"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gata_paper/CUTnTag_analysis/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
        -b)
            input_dir="$2"
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
	echo -e "Please provide the path to /Gata_paper/CUTnTag_analysis/ with -p"
  exit 1
fi

if [[ $input_dir == '' ]]
then
	echo -e "Please provide the path to a directory containing output from master1_CUTnTag_processing.sh with -b."
  exit 1
fi


input_dir=$(realpath $input_dir)'/'
work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'

# Calculate and plot correlation between CUT&Tag samples
echo -e "Calculate pearson correlation for CUT&Tag samples"
${path}scripts/CUTnTag_pearson.sh ${input_dir}final_bam/ $work_dir
Rscript ${path}scripts/CUTnTag_pearson.R $work_dir

#Find GATA motifs within CUT&Tag peaks
echo -e "Finding GATA motifs with FIMO"
${path}scripts/Find_motifs.sh $path ${input_dir}merged_peaks/ $work_dir

#Plot GATA CUT&Tag expression within peaks
echo -e "Quantify GATA signal within peaks"
Rscript ${path}scripts/quant_motifs.R $work_dir $input_dir

#Perform motif enrichment analysis using AME
echo -e "Enriching for JASPAR motifs using AME"
Rscript ${path}scripts/ame_file_prep.R $work_dir $input_dir
${path}scripts/run_ame.sh $path $work_dir
Rscript ${path}scripts/ame_fig.R $work_dir

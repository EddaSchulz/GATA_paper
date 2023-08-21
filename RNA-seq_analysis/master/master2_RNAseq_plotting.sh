#!/bin/bash
# Master script for RNAseq visualization (Published in vivo data and TX timecourse)
fastq_dir=''
work_dir=$(pwd)'/'
path=''

help() {
   echo "Plots heatmaps and lineplots for RNAseq data (in vivo and TX timecourse)."
   echo
   echo "Syntax: ./master2_RNAseq_plotting.sh [d|p|h]"
   echo "options:"
   echo "d     Provides working directory (Standard is current directory)."
   echo "h     Prints this help."
   echo "p     Provide path to /Gata_paper/RNA-seq_analysis/. [mandatory]"
   echo
}

parse_args() {
    case "$1" in
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


work_dir=$(realpath $work_dir)'/'
path=$(realpath $path)'/'


# Plots heatmap for TX timecourse data and GATA factor lineplot
echo -e "Plots TX timecourse data"
Rscript ${path}scripts/plot_TX.R $path $work_dir

# Plots heatmap for activators in vivo
echo -e "Plot activator preimplantation heatmap"
Rscript ${path}scripts/plot_activators_preimplant.R $path $work_dir

# Plots heatmap for GATA factors in vivo
echo -e "Plot GATA factor embryo heatmap"
Rscript ${path}scripts/gata_embryo_heatmap.R $path $work_dir

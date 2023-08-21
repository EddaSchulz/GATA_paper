#!/bin/bash
# Builds genome for bowtie2

path=$1
work_dir=$2
genome=$3 # Standard mm10 and B6/Cast masked genome is found at /NGS_alignment/files/

mkdir -p  ${work_dir}genome

bowtie2-build $genome ${work_dir}genome/mm10

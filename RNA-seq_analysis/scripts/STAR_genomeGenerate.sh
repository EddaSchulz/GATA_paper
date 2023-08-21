#!/bin/bash
# Builds genome for STAR

path=$1
work_dir=$2
genome=$3 # B6/Cast masked genome
gtf=$4

mkdir -p  ${work_dir}genome

STAR --runMode genomeGenerate --genomeDir ${work_dir}genome --genomeFastaFiles $genome --sjdbGTFfile $gtf

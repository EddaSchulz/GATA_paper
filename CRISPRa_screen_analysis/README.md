# CRISPRa screen analysis

## Description
This folder contains code that was used for the CRISPRa Screen analysis in Ravid-Lustig et al., Nature Cell Biology 2023

## Software dependencies and operating systems
In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- MAGeCK (v0.5.9.3) Model-based analysis of CRISPR screens available at "https://sourceforge.net/p/mageck/wiki/Home/"

R Scripts can be run using R (v4.2.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- gridExtra (v2.3)
- tidyverse (v2.0.0)

## Reproduce analysis
The raw sequencing data of the CRISPRa Screen needs to be retrieved from GEO and renamed according to "/CRISPRa_screen_analysis/files/rename_GSM.txt" (GSE194018). For paired samples, only read 1 should be used. The mm10 reference genome needs to be downloaded from ncbi and stored as "/CRISPRa_screen_analysis/files/mm10.fa" (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/).
Data analysis and plotting can be performed using the master scripts in the "Gata_paper/CRISPRa_screen_analysis/master/" directory. For both master scripts it is necessary to provide a path to "Gata_paper/CRISPRa_screen_analysis/" (with -p). They don't depend on each other.

- (i)   "master1_MAGeCK_analysis.sh": Conducts MAGeCK analysis for the CRISPRa Screen. Requires directory with FASTQ files from GEO (with -f).
- (ii)   "master2_plot_screen_data.sh: Generates plots for the CRISPRa screen, including quality control. Input files are already provided in "/Gata_paper/CRISPRa_screen_analysis/files/".

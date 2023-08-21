# RNA-seq processing and visualization

## Description
Processing of XX and XO TX1072 RNA-seq timecourse data (Pacini et al., 2021, Ravid-Lustig et al., 2023).
Visualization of the timecourse data and of in vivo RNAseq/pseudobulk scRNAseq (Deng et al., 2014, Zhang et al., 2017, Argelaguet et al., 2019). Deng 2014 and Zhang 2017 data was processed according to "https://github.com/EddaSchulz/Xert_paper/tree/main/NGS_alignment". Argelaguet 2019 data was retrieved from https://github.com/rargelaguet/scnmt_gastrulation.


## Software dependencies and operating systems

In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- STAR (v2.7.5a) aligner from "https://github.com/alexdobin/STAR"

R Scripts can be run using R (v4.2.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- gridExtra (v2.3)
- Rsubread (v2.12.3)
- Rgb (v1.7.5)
- RColorBrewer (v1.1-3)
- tidyverse (v2.0.0)


## Reproduce analysis
The alignment and data processing can be performed using the master scripts stored in the "RNA-seq_analysis/master/" directory. For both master scripts, it is necessary to provide a path to "Gata_paper/RNA-seq_analysis/" (with -p). " directory. For both master scripts, it is necessary to provide a path to "Gata_paper/CUTnTag_analysis/" (with -p). The B6/Cast-masked mm10 genome has to be prepared according to https://github.com/EddaSchulz/Pacini_paper and stored as "/RNA-seq_analysis/files/N_masked_B6_Cast.fa". The RDS file for scRNAseq data from Argelaguet et al., 2019 has to be downloaded from "https://github.com/rargelaguet/scnmt_gastrulation" and stored as "/files/SingleCellExperiment.rds"
- (i)   "master1_RNAseq_processing.sh": Builds masked mm10 genome, aligns and processes TX1072 timecourse RNA-seq data. The relevant FASTQ files have to be retrieved GEO and renamed according to "/files/rename_GSM.txt". Provide directory containing the FASTQ files (with -f).
- (ii)  "master2_RNAseq_plotting.sh": Performs visualization for RNA-seq timecourse and pseudobulk scRNAseq data (Deng et al., 2014, Zhang et al., 2017, Argelaguet et al., 2019).

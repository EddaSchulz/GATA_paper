# CUT&Tag alignment and analysis

## Description
CUT&Tag targeting GATA2/3 in TS cells and GATA4/6 in XEN cells (XX lines). Additionally, the data contains H3K27ac and H3K27me3 for both lines. The script performs alignment and analysis, including motif enrichment.

## Software dependencies and operating systems
This folder stores all the scripts used for read alignment and processing of the various CUT&Tag data.

In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- Bowtie2 (v2.3.5.1) aligner from "https://github.com/BenLangmead/bowtie2"
- Deeptools2 (v3.4.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- MACS2 (v2.1.2) PYTHON script from "https://github.com/macs3-project/MACS"
- Picard tools (v2.18.25) collection of JAVA scripts from "https://broadinstitute.github.io/picard/"
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- TrimGalore (v0.6.4) PERL script from "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"
- MEME Suite (v5.5.1) collection of tools for motif discovery and enrichment from "https://meme-suite.org/"

R Scripts can be run using R (v4.2.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- gridExtra (v2.3)
- tidyverse (v2.0.0)
- ggrepel (v0.9.3)
- Rsubread (v2.12.3)


## Reproduce analysis
The alignment and data processing can be performed using the master scripts stored in the "CUTnTag_analysis/master/" directory. For both master scripts, it is necessary to provide a path to "Gata_paper/CUTnTag_analysis/" (with -p). The mm10 reference genome needs to be downloaded from ncbi and stored as "/CUTnTag_analysis/files/mm10.fa" (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/). The B6/Cast-masked mm10 genome has to be prepared according to https://github.com/EddaSchulz/Pacini_paper and stored as "/CUTnTag_analysis/files/N_masked_B6_Cast.fa". 
- (i)   "master1_CUTnTag_processing.sh": Builds masked mm10 genome, aligns and processes CUT&Tag data, calls peaks and creates coverage tracks. The respective FASTQ files have to be downloaded from GEO and renamed according to "/files/rename_GSM.txt". Provide directory containing the FASTQ files (with -f).
- (ii)  "master2_CUTnTag_analysis.sh": Performs motif enrichment and correlation analysis of CUT&Tag data. Depends on "master1_CUTnTag_processing.sh" (directory containing output should be provided with -b).


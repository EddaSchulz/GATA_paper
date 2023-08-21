# Published ChIP-seq and ATAC-seq processing

## Description
Processing of 8-cell ATAC-seq data (Wu et al. 2016) and GATA6-FLAG ChIP-seq data in male ESCs with Gata6 overexpression (Wamaitha 2015).

## Software dependencies and operating systems

In order to perform these analyses, the following software has to be installed and available on the command line ($PATH):
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- Bowtie2 (v2.3.5.1) aligner from "https://github.com/BenLangmead/bowtie2"
- Deeptools2 (v3.4.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- MACS2 (v2.1.2) PYTHON script from "https://github.com/macs3-project/MACS"
- Picard tools (v2.18.25) collection of JAVA scripts from "https://broadinstitute.github.io/picard/"
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- TrimGalore (v0.6.4) PERL script from "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"


## Reproduce analysis
The alignment and data processing can be performed using the master scripts stored in the "ATAC-seq_Chip-seq_analysis/master/" directory. For both master scripts, it is necessary to provide a path to "Gata_paper/ATAC-seq_Chip-seq_analysis/" (with -p). The mm10 reference genome needs to be downloaded from ncbi and stored as "/CUTnTag_analysis/files/mm10.fa" (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/). For both master scripts, the relevant FASTQ files have to be retrieved GEO and renamed according to "/files/rename_GSM.txt". Provide directory containing the FASTQ files (with -f).
- (i)   "master1_ChIP_processing.sh": Builds mm10 genome, aligns and processes GATA6 ChIP-seq data from Wamaitha et al., 2015.
- (ii)  "master2_ATAC_Processing.sh": Builds mm10 genome, aligns and processes 8-cell embryo ATAC-seq data from Wu et al., 2016.

This repository contains all the code for processing and analysis of benzonase cross-linking-and-immunoprecipitation (bCLIP) for Integrator complex proteins in mouse embryonic stem cells (mESC) and RNA-sequencing (RNA-seq) data for siRNA-knockdowns of Integrator components in mESC and human HEK293 cell lines, generated in the lab of prof. Dr. Stefanie Jonas by Moes Murielle.
In addition, several public datasets were downloaded and analyzed alongside:

- enhanced cross-linking-and-immunoprecipitation (eCLIP) data for Ints11 in human Hela cells from Jasmine Barra et al (2020), https://www.science.org/doi/10.1126/sciadv.aaz9072
- bCLIP data from Matyas Flemr et al (2023), https://rnajournal.cshlp.org/content/29/8/1140.long
- PRO-seq data for Ints11 degron system in mESC, RNA-seq for Ints11 KD and Ints11 degron in mESC, TT-seq for Ints11 KD and Ints11 degron in mESC, Chip-seq for Ints11 in mESC from Chad B. Stein (2022), https://www.sciencedirect.com/science/article/pii/S1097276522009613?via%3Dihub

Briefly, the computations were done in the following steps:

**1)** Running jupyter notebook "bCLIP.ipynb".
It contains the code to download necesary public data, organize the input files into folders, create custom annotation .gtf files (union of GENCODE + RNAcentral), and prepare the configuration for the snakemake workflow.

**2)** Section "Run main WF" in the notebook "bCLIP.ipynb" creates the bash commands to run the respective snakemake workflows for bCLIP, eCLIP, and ChipSeq data (Snakefile_bCLIP), PRO-seq data (Snakefile_PROseq), RNA-seq and TT-seq data (Snakefile_RNAseq). The workflow files are located in the subfolder "./bclip_workflow".
Workflows should be run from within the subdirectory "./bclip_workflow".

**3)** The jupyter notebook "QC_plots.ipynb" produces tables and figures for the in-depth quality control analysis of bCLIP/eCLIP data. 

**4)** The jupyter notebook "QC_plots_detailed_mapping.ipynb" produces the figures to analyze the CIGAR and MD features from .bam files (soft-clip positions and locations of single-nucleotide polymorphisms within reads).

**5)** The jupyter notebook "Analysis_bCLIP.ipynb" contains all the code to generate final figures and tables for the manuscript.

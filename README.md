# DiscoVir: Reprodutable Viral Discovery Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Conda-Env](https://img.shields.io/badge/conda-env-green.svg)](workflow/envs/DiscoVir.yaml)
[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://opensource.org/licenses/GNU)

<p align="center">
  <img src="resources/logo/DiscoVir_Logo.png" width="300" alt="DiscoVir Logo">
</p>

---

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Author](#author)

## Introduction

DiscoVir is a comprehensive and scalable Snakemake workflow for the detection of novel viruses within HTS data. The pipeline is designed to work with both short-read (Illumina) and long-read (Nanopore/PacBio) metagenomic and transcriptomic data.

## Requirements

### Conda Installation
The pipeline uses Conda to manage its dependencies. We recommend the presence of any version running within the machine. If no version is available, follown the conda installation within miniforge, [conda](https://forge.ird.fr/transvihmi/nfernandez/install_conda_with_miniforge).

### Databases and Taxonomic files
Due the large size of databases used, and its common presence within HPC clusters, databases path must be prior present to run the pipeline and its download does not make part of this pipeline. In case of doubts, open an issue or get in touch within your HPC support team to proper intallation.
 - [Nr5](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/)
 - [prot.accession2taxid.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/)
 - [tamdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
 - [Kraken2](https://benlangmead.github.io/aws-indexes/k2)

## Installation
The DiscoVir.sh script is responsible to manage the run, automatically downloading and creating the necessary Conda environments, as well as manage conda versions for the pipeline to work.

1.  **Clone the repository:**

```bash
 git clone https://github.com/user/discovir.git
 cd discovir
```

2. **Verify installation:**

```shell
 bash DiscoVir.sh --help 
```

The following must appear within your screen.

```yaml
 DiscoVir: Viral Metagenomics & 'Dark Matter' Discovery

 Author: MSc. Matheus Cosentino 
 Version: 12.2025

 Usage:
  bash DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 Required Arguments:
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 Optional Arguments:
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)

 Database Overrides (Use external DBs):
   --diamond_db <FILE>  External Diamond database (.dmnd).
   --kraken2 <DIR>      External Kraken2 database directory (Must contain hash.k2d, opts.k2d, taxo.k2d)
   --taxdump <DIR>      Directory containing nodes.dmp and names.dmp

 Module Toggles (Enable/Disable Analysis):
   --assembly           Enable De Novo Assembly (Default: False)
   --kraken2            Enable Kraken2 Taxonomy (Default: False)
   --diamond            Enable Diamond Taxonomy (Default: False)
   --darkmatter         Enable Palm Annot / Dark Matter (Default: False)
   --remove-download    Revome downloaded SRA files (Default: Keep)
   --skip-reads-kraken  Skip Reads Kraken2 ID (Default: Run)
   --skip-reads-diamond Skip Reads Diamond ID (Default: Run)
 
 Flags:
   -h, --help           Show this help message
   -v, --version        Show version
```


## Usage
### Overall Pipeline Use
The DiscoVir standard configuration is optimized for Illumina data and to run within an HPC cluster managed by [slurm](https://slurm.schedmd.com/documentat). The script is intrinsically capable of automatically differentiating between paired-end and single-end data formats as well as to diferentiate within fasta format. - Kraken2 & Diamond Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --kraken2 <DIR> --diamond_db <FILE> --taxdump <DIR> 
```

- Kraken2 Only Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR>  --kraken2 <DIR> --skip-reads-diamond
```

- Diamond Only Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --diamond_db <FILE> --taxdump <DIR> --skip-reads-kraken
```

### DiscoVir Core Pipelline
After a first round of analysis that the sequencing data demonstrates an overall satisfatory quality, the following command can be used to a further viral discovery pipelline. This process consists in a contig assembly within spades viralrna aalgorythim (or any other tag or assembler of your preference among options in the config.yaml file), followed by taxonomic identification within diamond. Contigs without diamond hits are filtered by size and Aminoacid generated by ORFs within distinct genetic codes pass though the [palm_annot](https://github.com/rcedgar/palm_annot.git) RdRp scan scripts. Total Virus contigs iddentified by viral family are reported within a Rmarkdom script.

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --diamond --duskmatter 
```

## Features
### Working with Public Data

One of the advantages within the DiscoVir pipeline, is the possibility to easily incorporate public available data within the discovery pipeline, allowing an easy process in discovering novel viruses within public [SRA](https://www.ncbi.nlm.nih.gov/sra/).

For it, an line separated file (accession.txt) file must be provided, where each line corresponds to a Run iddentifier within SRA.

```yaml
SRR25008973
SRR25008974
SRR25008975
SRR25008976
SRR25008977
```

Once librarys of interest are givenn, run the following command to download the files.

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --sra <FILE> 
```

**Caution**:
After the whole workflow is finished, data will be keepted. In case of interested in deleting the data after the run, use the following command:

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --sra <FILE> --remove-download 
```

### Oxford Nanopore Sequencing Data

The proposed workflow is capable of dealing with long-read data, such as the HTS files generated by ONP. However, editions to configuration file must be made due the variation within the type of Reads.

- Adjust the overall quality and size filtering criteria, modify the relevant parameters for fastp inside the config/config.yaml file.

```yaml
# QC
fastp:
  length_required: 
    #- 15 #illlumina adapters expected
    - 500 #nanopore conservative size for low reads
  qualified_quality_phred: 
    #- 30 #illumina paired end
    - 10 #long read
```

After edition, run the commands expected for read analysis according to your interest.

- Kraken2 & Diamond Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --kraken2 <DIR> --diamond_db <FILE> --taxdump <DIR> 
```

- Kraken2 Only Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR>  --kraken2 <DIR> --skip-reads-diamond
```

- Diamond Only Reads Iddentfication
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --diamond_db <FILE> --taxdump <DIR> --skip-reads-kraken
```

### ONP De Novo Assembly

In scenarios where de novo assembly is necessary for ONP data, the following alterations within the config must be done:
 - Activate the tools of interest. Medaka contigs correction is recommended, but models of correction varys according to flowncell and sequencer used.

```yaml
tool:
  denovo: #select your favorite de novo tools
  # options: 
   # - 'spades'    # more resourcers needed, but algorythis to better deal with RNA seq data of Illumina short reads(default)
   # - 'megahit'    # indicated for llumina short reads. Faster, but more error proner
    - 'flye'       #indicated for nanopore sequencing long reads (>1000 bp per reads)
   # - 'raven'      #indicated for nanopore sequencing. Faster, but more error proner
    - 'medaka_flye'   #correction to ONP flye assembly
   # - 'medaka_raven'  #correction to ONP raven assembly
```

The same process can be made to use the raven assemby software and medaka correction, also available. Nevertheless, is worth mentioning that overall contig assembly within ONP reads overall generate very low number of contigs.


# Contributing
Contributions are welcome! If you find a bug or have a suggestion for improvement, please open an issue!

# License
This project is licensed under the GNU General Public License v3.0.

# Author
MSc. Matheus Cosentino

---
<div align="center">
  <img src="resources/logo/LDDV.png" width="200" alt="LDDV Logo" style="display: inline-block; margin-left: 10px;">
  <img src="resources/logo/index.png" width="200" alt="Transvhmi Logo" style="display: inline-block; margin-left: 10px;">
</div>
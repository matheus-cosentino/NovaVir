# DiscoVir: Reproducible Viral Discovery Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Conda-Env](https://img.shields.io/badge/conda-env-green.svg)](workflow/envs/DiscoVir.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0)

<p align="center">
  <img src="resources/logo/DiscoVir_Logo.png" width="300" alt="DiscoVir Logo">
</p>

---

## Table of Contents

- [Introduction](#introduction)
- [Workflow Overview](#workflow-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Expected Outputs](#expected-outputs)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Author](#author)

## Introduction

DiscoVir is a comprehensive and scalable Snakemake workflow for the detection of novel viruses within High-Throughput Sequencing (HTS) data. The pipeline is designed to work with both short-read (Illumina) and long-read (Nanopore/PacBio) metagenomic and transcriptomic data.

## Workflow Overview

```mermaid
graph TD;
    A[Raw Reads fastq.gz] --> B(fastp Quality Control);
    B -->|Filtered Reads| C{Analysis Path};
    C -->|Read-Level Taxonomy| D[Kraken2 & Diamond];
    C -->|De Novo Assembly| E[SPAdes / Megahit / Flye];
    
    E -->|Contigs| F(Diamond Annotation vs NR);
    E -->|Contigs| G(Kraken2 Annotation);
    E -->|Contigs| H(Dark Matter Module);
    
    H -->|ORF Prediction & HMM| I(PalmAnnot2 & RVDB Screening);
    I -->|Taxonomy Summarization| J[Final RdRp Hits & Annotations];

    D --> K[Krona & Report Dashboards];
    F --> K;
    J --> K;
```

## Requirements

### Recommended Hardware Specs
- **Memory (RAM)**: Minimum 32GB (for metaSPAdes assembly). 64GB+ recommended for large meta-transcriptomic datasets.
- **CPU**: 8-16+ cores recommended.
- **Disk Space**: At least 100GB+ available space to download external databases safely on local or HPC clusters.

### Conda Installation
The pipeline uses Conda to manage its dependencies. We recommend having any version of Conda running on your machine. If you don't have Conda installed, you can follow the installation guide for Miniforge at [conda-forge](https://github.com/conda-forge/miniforge).

### Databases
Due to the large size of the databases used and their common presence within HPC clusters, the paths to the databases must be provided before running the pipeline. Their download is not part of this pipeline. If you have any doubts, please open an issue or get in touch with your HPC support team for proper installation.
 - [nr](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/) (For Diamond)
 - [Kraken2](https://benlangmead.github.io/aws-indexes/k2)

## Installation
The `DiscoVir.sh` script is responsible for managing the run, automatically downloading and creating the necessary Conda environments, as well as managing Conda versions for the pipeline to work.

1.  **Clone the repository:**

```bash
 git clone https://github.com/matheus-cosentino/discovir.git
 cd discovir
```

2. **Verify installation:**

```shell
 bash DiscoVir.sh --help 
```

The following must appear on your screen:

 ```yaml
 DiscoVir: Viral Metagenomics & 'Dark Matter' Discovery

 Author: MSc. Matheus Cosentino 
 Version: 03.2026

 Usage:
  bash DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 Required Arguments:
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 Optional Arguments:
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)
   --temp-dir <DIR>     Temporary directory (default: /tmp)

 Database Overrides (Use external DBs):
   --diamond_db <FILE>  External Diamond database (.dmnd).
   --kraken2 <DIR>      External Kraken2 database directory (Must contain hash.k2d, opts.k2d, taxo.k2d)

 Module Toggles (Enable/Disable Analysis):
   
- **`--reads-kraken`**: Enables taxonomic classification of raw reads using Kraken2.
- **`--reads-diamond`**: Enables taxonomic classification of raw reads using DIAMOND.
- **`--assembly`**: Enables *de novo* assembly of reads into contigs.
- **`--kraken2`**: Enables taxonomic classification of assembled contigs using Kraken2.
- **`--diamond`**: Enables taxonomic classification of assembled contigs using DIAMOND.
- **`--darkmatter`**: Enables the "dark matter" module to search for RdRp signatures & RVDB evaluation.
- **`--remove-download`**: Removes SRA files after download.

 
 Flags:
   -h, --help           Show this help message
   -v, --version        Show version
```

## Usage

### Core Concepts for First-Time Users

DiscoVir's design is highly modular. The use of specific **flags** triggers different modules with distinct purposes and computational requirements:

- **Initial Evaluation (`--reads-kraken`, `--reads-diamond`)**: These modules operate directly on raw reads. They serve as a rapid, initial evaluation of your pipeline and sample composition, allowing you to gauge the presence of viral content without heavy computation.
- **Deep Viral Discovery (`--kraken2`, `--diamond`, `--darkmatter`)**: These modules heavily depend on the prior *de novo* assembly of your reads (`--assembly` is usually implicitly run or required). Because they operate on assembled contigs and search for novel sequences, they focus on more computational power and time, yielding highly specific viral insights.
- **Resource Management (`--profile`)**: **The profile flag is key** to running DiscoVir successfully. It identifies the type of computational resources (CPU, Memory, Time) required for each specific rule. Properly setting your `--profile` (e.g., local execution vs. HPC Slurm) is essential for a smooth first-time experience.

### Example Workflows

- **Initial Evaluation: Kraken2 & Diamond Reads Identification**
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --profile <STR> --kraken2 <DIR> --diamond_db <FILE> --reads-kraken --reads-diamond
```

- **Initial Evaluation: Kraken2 Only Reads Identification**
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --profile <STR> --kraken2 <DIR> --reads-kraken
```

- **Initial Evaluation: Diamond Only Reads Identification**
```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --profile <STR> --diamond_db <FILE> --reads-diamond
```

### DiscoVir Core Pipeline (Assembly & Discovery)
After a rapid initial analysis demonstrates satisfactory overall quality, the following command can be used for the deep viral discovery pipeline. This process requires more computational power and consists of a contig assembly with the SPAdes `viralrna` algorithm (or any requested assembler in `config.yaml`), followed by taxonomic identification. 

Contigs without Diamond hits are filtered by size, and Amino Acids generated by ORFs within distinct genetic codes are passed through the [palm_annot](https://github.com/rcedgar/palm_annot.git) RdRp scan scripts along with rapid [RVDB HMM Summaries](https://rvdb-prot.pasteur.fr/). Viral contigs identified by viral family are reported within an R Markdown script, alltogether with palm_annot results. RVDB summaries are in development.

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --profile <STR> --diamond --darkmatter 
```

## Expected Outputs

Once successfully executed, DiscoVir cleanly organizes your results within your standard `--output` directories. Below is a generic overview of an assembled mapping workflow:

```text
results/
├── log/                           # Main execution logs
├── SRR10677983/
│   ├── fastp/                     # QC & Trimmed adapter reports
│   ├── spades/                    # Assembled Fasta contigs
│   ├── diamond_spades_kauto/      # Annotated Diamond outputs vs NR
│   ├── kraken2_spades_kauto/      # Fast Kraken2 K-mer validations
│   ├── rvdb_spades_kauto/         # RVDB Summaries containing precise NCBI translation
│   ├── darkmatter_spades_kauto/   # ORFs & Filtered RdRp viral candidates!
│   └── krona_spades_kauto/        # Interactive HTML plots of relative abundances
└── multiqc_report.html            # Aggregated sequencing QC metrics dashboard
```


## Features
### Working with Public Data

One of the advantages of the DiscoVir pipeline is the possibility to easily incorporate publicly available data into the discovery pipeline, allowing for an easy process of discovering novel viruses within the public [SRA](https://www.ncbi.nlm.nih.gov/sra/).

For this, a line-separated file (e.g., `accession.txt`) must be provided, where each line corresponds to a Run identifier within SRA.

```yaml
SRR25008973
SRR25008974
SRR25008975
SRR25008976
SRR25008977
```

Once the libraries of interest are given, run the following command to download the files.

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --sra <FILE> 
```

**Caution**:
After the whole workflow is finished, the data will be kept. If you are interested in deleting the data after the run, use the following command:

```shell
bash DiscoVir.sh --input <DIR> --output <DIR> --sra <FILE> --remove-download 
```

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
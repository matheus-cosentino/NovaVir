# DiscoVir: Scalable Viral Discovery Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Conda-Env](https://img.shields.io/badge/conda-env-green.svg)](workflow/envs/DiscoVir.yaml)
[![License: GNU](https://img.shields.io/badge/License-GNU-yellow.svg)](https://opensource.org/licenses/GNU)

<p align="center">
  <img src="resources/logo/DiscoVir_Logo.png" width="300" alt="DiscoVir Logo">
</p>

<p align="center">
  <strong>Uncovering viral diversity, evolutionary dynamics, and ecological patterns from paired and unpaired reads data using Snakemake and SLURM.</strong>
</p>

---

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output](#output)
- [Contributing](#contributing)
- [License](#license)
- [Author](#author)

## Introduction

DiscoVir is a comprehensive and scalable Snakemake workflow for the detection of known viruses and the discovery of novel members of the viral "dusk matter". It leverages the RNA-dependent RNA polymerase (RdRp) barcode for identification. The pipeline is designed to work with both short-read (Illumina) and long-read (Nanopore/PacBio) metagenomic and transcriptomic data.

## Features

DiscoVir is a modular and flexible pipeline. The main modules are:

-   **Data Input & Quality Control**:
    -   Supports raw reads (.fastq.gz) or pre-assembled contigs (.fasta).
    -   Automatically fetches samples from the Sequence Read Archive (SRA) if an accession list is provided.
    -   Performs trimming and quality filtering using `fastp`.
-   **De Novo Assembly**:
    -   Supports multiple assemblers for different sequencing technologies:
        -   **SPAdes**: Optimized for Illumina/RNA viral data.
        -   **Megahit**: A fast and memory-efficient assembler for large metagenomes.
        -   **Flye**: Specialized for long-read data (Nanopore/PacBio).
-   **Taxonomic Classification**:
    -   **Kraken2**: High-speed taxonomic assignment using k-mers.
    -   **Diamond**: Sensitive protein alignment against the NR database.
    -   **TaxonKit**: For lineage resolution of identified hits.
-   **"Dusk Matter" Discovery (Palm Annot)**:
    -   Targeted search for RNA-dependent RNA polymerase (RdRp) motifs using `PalmDB`. This module is crucial for identifying divergent RNA viruses that lack homology in standard databases.
-   **Reporting**:
    -   Generates a final HTML report summarizing QC, assembly statistics, and viral findings for every sample.

## Installation

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/user/discovir.git
    cd discovir
    ```

2.  **Install Conda and Mamba:**

    The pipeline uses Conda to manage its dependencies. We recommend using [Mamba](https://github.com/mamba-org/mamba) for faster environment creation.

    ```bash
    conda install -c conda-forge mamba
    ```

3.  **Create the main Conda environment:**

    The `DiscoVir.sh` script will automatically create the necessary Conda environments for the pipeline. However, you can create the main environment manually:

    ```bash
    mamba env create -f workflow/envs/DiscoVir.yaml
    ```

## Configuration

The pipeline is configured using the `config/config.yaml` file. Here you can define input/output paths, select modules to run, and set parameters for each tool.

**Main options:**

-   `input_dir`: Path to the directory containing your input files.
-   `output_dir`: Path to the directory where results will be saved.
-   `sra_accessions`: A file containing SRA accession numbers to download.
-   `modules`: A boolean flag to enable or disable specific modules (`assembly`, `kraken2`, `diamond`, `duskmatter`).
-   `tool`: Select the assembler to use (`spades`, `megahit`, `flye`).

## Usage

The pipeline is wrapped in a Bash script (`DiscoVir.sh`) that handles Conda environment creation, resource linking, and job submission. DiscoVir is optimized for SLURM environments. To run on a cluster, you can specify the number of parallel jobs:

```bash
bash DiscoVir.sh --input data/raw --output results_project --jobs 50
Arguments
--input: Directory containing raw reads (.fastq.gz) or contigs (.fasta).
--output: Directory where results and reports will be saved.
--sra <FILE>: Text file containing SRA Accession IDs for automated download.
--jobs <INT>: Number of parallel jobs (default: 15).
--kraken2 <DIR>: Path to an external Kraken2 database.
--diamond_db <FILE>: Path to an external Diamond .dmnd database.
```

Output
The output directory will be structured as follows:

```bash
<output_dir>/
├───0_raw_reads/
├───1_preprocessed_reads/
├───2_assembly/
├───3_kraken2/
├───4_diamond/
├───5_duskmatter/
└───report/

raw_reads/: Raw input reads (or symlinks to them).
filteres_reads/: Quality-controlled reads from fastp.
assembly/: Assembled contigs from SPAdes, Megahit, or Flye.
kraken2/: Taxonomic classification results from Kraken2.
diamond/: Protein alignment results from Diamond.
duskmatter/: RdRp annotation results.
report/: Final HTML report & Suplementary Files
```


Contributing
Contributions are welcome! If you find a bug or have a suggestion for improvement, please open an issue or submit a pull request.

License
This project is licensed under the GNU General Public License v3.0. See the LICENSE">LICENSE file for details.

Author
MSc. Matheus Cosentino


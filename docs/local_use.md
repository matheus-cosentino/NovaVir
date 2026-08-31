# NovaVir: First Time Local Usage

This document guides you through the initial setup required to run NovaVir on your local machine, including installing Conda, preparing testing databases, and executing the pipeline.

## 1. Getting Conda

NovaVir heavily relies on Conda for dependency management, ensuring reproducibility across environments. The pipeline uses Snakemake's native `--use-conda` flag to automatically create environments for each tool.

If you don't have Conda installed, we recommend installing **Miniforge** (which uses conda-forge as the default channel and is free for all usage):
1. Visit the [Miniforge GitHub repository](https://github.com/conda-forge/miniforge).
2. Download the installer for your operating system.
3. Follow the installation instructions provided on the repository.

Once installed, ensure `conda` is available in your terminal:
```bash
conda --version
```

## 2. Quick Start: Test Databases

If you are just testing the pipeline and don't want to download massive databases right away, you can use the *E. coli* K-12 proteome as a quick testing database for DIAMOND, and the RefSeq Viral database for Kraken2.

### Step 2.1: E. coli K-12 DIAMOND Database

```bash
# Create directory
mkdir -p db/diamond && cd db/diamond

# Baixa o proteoma da E. coli K-12
wget "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:83333)" -O ecoli_proteome.fasta 

# Formata o cabeçalho para manter apenas a ID da proteína
awk -F'|' '/^>/ {print ">"$2} !/^>/ {print $0}' ecoli_proteome.fasta > ecoli_formatted.fasta

# Constrói o banco do DIAMOND
diamond makedb --in ecoli_formatted.fasta -d ecoli_db

# Remove os arquivos FASTA temporários
rm -rf ecoli_proteome.fasta ecoli_formatted.fasta
cd ../../
```

### Step 2.2: Kraken2 RefSeq Viral Database

For Kraken2, you can use the lightweight RefSeq Viral database provided by the [Kraken2 AWS Indexes](https://benlangmead.github.io/aws-indexes/k2).

```bash
mkdir -p db/kraken2_viral && cd db/kraken2_viral

# Download the RefSeq Viral database (approx 0.5 GB)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz

# Extract the database
tar -zxvf k2_viral_20240112.tar.gz

# Clean up
rm k2_viral_20240112.tar.gz
cd ../../
```

---

## 3. Running NovaVir Locally

Once your databases are prepared, you can launch the pipeline on your local machine. NovaVir uses Snakemake profiles to manage resources. For local execution, we use the `profiles/local` profile.

### Step 3.1: Prepare Your Input Directory

Create a directory containing your raw sequencing reads (e.g., `.fastq.gz` files).

```bash
mkdir -p my_project/raw_reads
# Copy your fastq files into my_project/raw_reads
cp sample1_R1.fastq.gz my_project/raw_reads/
cp sample1_R2.fastq.gz my_project/raw_reads/
```

### Step 3.2: Execute the Pipeline

Run the `NovaVir.sh` script, pointing it to your local profile, your input directory, and the test databases you just created.

```bash
# Return to the NovaVir repository root directory
cd /path/to/novavir

# Run the pipeline locally
bash NovaVir.sh \
    --input ../my_project/raw_reads \
    --output ../my_project/results \
    --profile profiles/local \
    --diamond_db db/diamond/ecoli_db.dmnd \
    --kraken2 db/kraken2_viral \
    --diamond \
    --kraken2_reads \
    --darkmatter \
    --jobs 8
```

### Explanation of Arguments:

- `--input`: Path to the directory containing your `.fastq.gz` or `.fasta` files.
- `--output`: Path where the pipeline results will be saved.
- `--profile profiles/local`: Tells Snakemake to run locally and adhere to the resource limits defined in the local profile.
- `--diamond_db`: The exact path to the `ecoli_db.dmnd` database you compiled.
- `--kraken2`: The directory containing the extracted Kraken2 RefSeq Viral database.
- `--diamond`: Enables the assembly-based DIAMOND classification module.
- `--kraken2_reads`: Enables the reads-based Kraken2 classification module.
- `--darkmatter`: Enables the RdRp 'dark matter' discovery module.
- `--jobs 8`: Maximum number of parallel cores/threads Snakemake is allowed to use.

### Step 3.3: Reviewing the Results

When the pipeline finishes, your results will be structured within the specified `--output` directory. Key files include:

- `results/<sample>/darkmatter_spades_kauto/<sample>_Report_Diversity.html`: The interactive HTML report summarizing the viral diversity and RdRp candidates.
- `results/<sample>/diamond_spades_kauto/<sample>_spades_kauto_hits_with_lineage.tsv`: The full DIAMOND output augmented with complete NCBI lineages.
- `results/<sample>/krona_spades_kauto/`: Contains interactive Krona charts for navigating the taxonomic composition.

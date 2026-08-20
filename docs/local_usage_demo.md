# DiscoVir: Local Usage Demonstration

This document provides a step-by-step demonstration of how to run the DiscoVir pipeline locally. It covers everything from preparing the DIAMOND database with the required `>ptn_acc` header format to executing the pipeline on your local machine.

## 1. Preparing the DIAMOND Database

To perform accurate taxonomic classification with DIAMOND and MEGAN, the database must be correctly formatted and linked to the NCBI taxonomy. 

DIAMOND requires that the FASTA headers contain only the protein accession number (the `>ptn_acc` format).

### Step 1.1: Download the NR Database and Taxonomy Files

First, download the NCBI non-redundant (NR) protein database and the necessary taxonomy mapping files.

```bash
# Create a directory for the database
mkdir -p db/diamond && cd db/diamond

# Download the NR fasta file
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

# Download the taxonomy nodes and names
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip nodes.dmp names.dmp

# Download the accession to taxid mapping file
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

### Step 1.2: Format the FASTA Headers to `>ptn_acc`

NCBI FASTA headers typically contain descriptions along with the accession number (e.g., `>NP_047200.1 Replicase [Virus...]`). To meet the `>ptn_acc` format requirement, we must remove everything after the first space.

You can use `awk` to quickly format the FASTA headers:

```bash
# Keep only the accession number in the FASTA header
awk '{print $1}' nr > nr_formatted.fasta
```

The headers in `nr_formatted.fasta` will now look like:
```text
>NP_047200.1
>YP_009137151.1
```

### Step 1.3: Build the DIAMOND Database

Now, compile the DIAMOND database using the formatted FASTA file and the downloaded taxonomy files. This allows DIAMOND to automatically assign Taxonomic IDs to each hit.

```bash
diamond makedb \
    --in nr_formatted.fasta \
    --db nr \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp \
    --taxonnames names.dmp \
    --threads 8
```

This will output a compiled database file named `nr.dmnd`, ready to be used by DiscoVir.

---

## 2. Running DiscoVir Locally

Once your databases are prepared, you can launch the pipeline on your local machine. DiscoVir uses Snakemake profiles to manage resources. For local execution, we use the `profiles/local` profile.

### Step 2.1: Prepare Your Input Directory

Create a directory containing your raw sequencing reads (e.g., `.fastq.gz` files).

```bash
mkdir -p my_project/raw_reads
# Copy your fastq files into my_project/raw_reads
cp sample1_R1.fastq.gz my_project/raw_reads/
cp sample1_R2.fastq.gz my_project/raw_reads/
```

### Step 2.2: Execute the Pipeline

Run the `DiscoVir.sh` script, pointing it to your local profile, your input directory, and the external DIAMOND database you just created.

```bash
# Return to the DiscoVir repository root directory
cd /path/to/discovir

# Run the pipeline locally
bash DiscoVir.sh \
    --input ../my_project/raw_reads \
    --output ../my_project/results \
    --profile profiles/local \
    --diamond_db ../db/diamond/nr.dmnd \
    --diamond \
    --darkmatter \
    --jobs 8
```

### Explanation of Arguments:

- `--input`: Path to the directory containing your `.fastq.gz` or `.fasta` files.
- `--output`: Path where the pipeline results will be saved.
- `--profile profiles/local`: Tells Snakemake to run locally and adhere to the resource limits defined in the local profile.
- `--diamond_db`: The exact path to the `nr.dmnd` database you compiled.
- `--diamond`: Enables the assembly-based DIAMOND classification module.
- `--darkmatter`: Enables the RdRp 'dark matter' discovery module.
- `--jobs 8`: Maximum number of parallel cores/threads Snakemake is allowed to use.

### Step 2.3: Reviewing the Results

When the pipeline finishes, your results will be structured within the specified `--output` directory. Key files include:

- `results/<sample>/darkmatter_spades_kauto/<sample>_Report_Diversity.html`: The interactive HTML report summarizing the viral diversity and RdRp candidates.
- `results/<sample>/diamond_spades_kauto/<sample>_spades_kauto_hits_with_lineage.tsv`: The full DIAMOND output augmented with complete NCBI lineages.
- `results/<sample>/krona_spades_kauto/`: Contains interactive Krona charts for navigating the taxonomic composition.

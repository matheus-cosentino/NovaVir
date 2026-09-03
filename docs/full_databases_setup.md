# Full Databases (For Production)

For deep taxonomic classification and real viral discovery, you will need the full databases. Due to their large size, they are not bundled with DeepVir and must be downloaded and prepared beforehand. If you are on an HPC cluster, check with your administrators as these databases might already be available.

You will need:
1. **Kraken2 Database**: Can be downloaded from the [Kraken2 AWS Indexes](https://benlangmead.github.io/aws-indexes/k2). We recommend the "Standard" or "PlusPFP" database depending on your available memory (Standard requires ~50GB RAM, PlusPFP requires ~80GB RAM).
2. **DIAMOND Database (NR)**: You need the NCBI non-redundant (NR) protein database for deep taxonomic classification.

## Formatting the Full DIAMOND Database

To perform accurate taxonomic classification with DIAMOND and MEGAN, the NR database must be correctly formatted and linked to the NCBI taxonomy. DIAMOND requires that the FASTA headers contain **only the protein accession number** (the `>ptn_acc` format).

### Step 1: Download the NR Database and Taxonomy Files

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

### Step 2: Format the FASTA Headers to `>ptn_acc`

NCBI FASTA headers typically contain descriptions along with the accession number (e.g., `>NP_047200.1 Replicase [Virus...]`). To meet the `>ptn_acc` format requirement, we must remove everything after the first space.

You can use `awk` to quickly format the FASTA headers:

```bash
# Keep only the accession number in the FASTA header
awk '{print $1}' nr > nr_formatted.fasta
```

### Step 3: Build the DIAMOND Database

Compile the DIAMOND database using the formatted FASTA file and the downloaded taxonomy files. This allows DIAMOND to automatically assign Taxonomic IDs to each hit.

```bash
diamond makedb \
    --in nr_formatted.fasta \
    --db nr \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp \
    --taxonnames names.dmp \
    --threads 8
```

This will output a compiled database file named `nr.dmnd`, which you will pass to DeepVir using the `--diamond_db` flag.

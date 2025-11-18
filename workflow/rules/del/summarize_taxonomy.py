# workflow/scripts/summarize_taxonomy.py

import pandas as pd
import gzip
import sys

# --- Snakemake I/O ---
hits_file = snakemake.input.hits
nodes_file = snakemake.input.nodes
names_file = snakemake.input.names
accession_map_file = snakemake.input.accession_map
output_file = snakemake.output.summary
desired_rank = snakemake.params.rank

print(f"Processing DIAMOND hits file: {hits_file}", file=sys.stderr)
# --- 1. Get the list of protein accessions we need to look up ---
# Use a set for fast, unique storage
needed_accessions = set(pd.read_csv(hits_file, sep='\t', header=None, usecols=[1])[1].unique())
print(f"Found {len(needed_accessions)} unique protein accessions to map.", file=sys.stderr)

# --- 2. Build a mapping dict for only the accessions we need ---
print(f"Building accession-to-taxid map from {accession_map_file}...", file=sys.stderr)
accession_to_taxid = {}
# The accession file has NO version number, while DIAMOND's might. We check both.
with gzip.open(accession_map_file, 'rt') as f:
    # Skip header
    next(f)
    for line in f:
        accession_version, accession, tax_id, _ = line.strip().split('\t')
        if accession in needed_accessions or accession_version in needed_accessions:
            accession_to_taxid[accession] = int(tax_id)
            accession_to_taxid[accession_version] = int(tax_id)

print(f"Successfully mapped {len(accession_to_taxid)} accessions.", file=sys.stderr)

# --- 3. Load NCBI Taxonomy Data ---
print("Loading NCBI taxonomy nodes and names...", file=sys.stderr)
names = {}
with open(names_file, 'r') as f:
    for line in f:
        if 'scientific name' in line:
            parts = line.split('\t|\t')
            tax_id, name = parts[0], parts[1]
            names[int(tax_id)] = name

parents = {}
ranks = {}
with open(nodes_file, 'r') as f:
    for line in f:
        parts = line.split('\t|\t')
        tax_id, parent_id, rank = int(parts[0]), int(parts[1]), parts[2]
        parents[tax_id] = parent_id
        ranks[tax_id] = rank

def get_lineage_at_rank(tax_id, desired_rank):
    """Walks up the taxonomy tree to find a specific rank for a given tax_id."""
    current_id = tax_id
    for _ in range(30): # Safety break to prevent infinite loops
        if current_id not in parents or current_id not in ranks:
            return "Unclassified"
        if ranks[current_id] == desired_rank:
            return names.get(current_id, "Unknown Name")
        if current_id == 1: # Reached the root
            break
        current_id = parents[current_id]
    return "Unclassified"

# --- 4. Process DIAMOND Hits and Map to Taxa ---
print("Mapping hits to taxonomy...", file=sys.stderr)
df_hits = pd.read_csv(hits_file, sep='\t', header=None, usecols=[0, 1], names=['qseqid', 'sseqid'])
df_hits.drop_duplicates(inplace=True)

# Map protein accession to taxid
df_hits['staxids'] = df_hits['sseqid'].map(accession_to_taxid)
df_hits.dropna(subset=['staxids'], inplace=True)
df_hits['staxids'] = df_hits['staxids'].astype(int)

# Map taxid to the desired rank
print(f"Finding rank: {desired_rank}...", file=sys.stderr)
df_hits['Taxon'] = df_hits['staxids'].apply(lambda x: get_lineage_at_rank(x, desired_rank))

# --- 5. Summarize and Save ---
print("Counting taxa and saving summary...", file=sys.stderr)
summary_df = df_hits['Taxon'].value_counts().reset_index()
summary_df.columns = ['Taxon', 'Count']
summary_df.to_csv(output_file, sep='\t', index=False)

print("Done.", file=sys.stderr)

import json
import pandas as pd
import os
from datetime import datetime

# --- Snakemake I/O ---
output_md = snakemake.output.report_md
logo_path = snakemake.input.logo
report_title = snakemake.params.title
report_author = snakemake.params.author
fastp_json_files = snakemake.input.fastp_reports
spades_contigs_files = snakemake.input.spades_contigs
read_hits_files = snakemake.input.read_hits
contig_hits_files = snakemake.input.contig_hits
contig_summary_files = snakemake.input.contig_summaries

# --- Helper Function ---
def count_fasta_records(fasta_file):
    """Counts the number of records (sequences) in a FASTA file."""
    try:
        with open(fasta_file, 'r') as f:
            return sum(1 for line in f if line.startswith('>'))
    except FileNotFoundError:
        return 0

# --- Main Logic ---
with open(output_md, 'w') as report:
    # --- 1. Write the YAML Metadata Block for the PDF Title Page ---
    report.write("---\n")
    report.write(f"title: \"{report_title}\"\n")
    report.write(f"author: \"{report_author}\"\n")
    report.write(f"date: \"{datetime.now().strftime('%B %d, %Y')}\"\n")
    report.write(f"logo: {logo_path}\n")
    report.write("---\n\n")

    # --- 2. Write the rest of the report content ---
    report.write("# Workflow Summary\n\n")
    
    # Process each sample, assuming file lists are in the same order
    for i, fastp_json in enumerate(fastp_json_files):
        sample_name = os.path.basename(fastp_json).split('.')[0]
        report.write(f"## Sample: {sample_name}\n\n")
        
        # Read Processing Stats
        report.write("### Read Processing\n\n")
        with open(fastp_json, 'r') as f:
            stats = json.load(f)
            total_reads = stats['summary']['before_filtering']['total_reads']
            filtered_reads = stats['summary']['after_filtering']['total_reads']
            report.write(f"- **Total Raw Reads:** {total_reads:,}\n")
            report.write(f"- **Reads Kept After Filtering:** {filtered_reads:,}\n\n")

        # Assembly Stats
        report.write("### Assembly\n\n")
        num_contigs = count_fasta_records(spades_contigs_files[i])
        report.write(f"- **Contigs Generated:** {num_contigs:,}\n\n")
        
        # Taxonomic Identification Stats
        report.write("### Taxonomic Identification\n\n")
        try:
            read_hits_df = pd.read_csv(read_hits_files[i], sep='\t', header=None)
            identified_reads = read_hits_df[0].nunique()
            report.write(f"- **Reads Taxonomically Identified:** {identified_reads:,}\n")
        except (FileNotFoundError, pd.errors.EmptyDataError):
            report.write("- **Reads Taxonomically Identified:** 0\n")

        try:
            contig_hits_df = pd.read_csv(contig_hits_files[i], sep='\t', header=None)
            identified_contigs = contig_hits_df[0].nunique()
            report.write(f"- **Contigs Taxonomically Identified:** {identified_contigs:,}\n\n")
        except (FileNotFoundError, pd.errors.EmptyDataError):
            report.write("- **Contigs Taxonomically Identified:** 0\n\n")

        # Viral Diversity Summary from Contigs
        report.write("### Viral Diversity in Contigs\n\n")
        try:
            contig_summary_df = pd.read_csv(contig_summary_files[i], sep='\t')
            # Filter for common viral family names
            viral_df = contig_summary_df[contig_summary_df.iloc[:,0].str.contains('viridae', case=False)]
            if not viral_df.empty:
                report.write(viral_df.to_markdown(index=False))
            else:
                report.write("No viral families identified in contigs.")
            report.write("\n\n")
        except (FileNotFoundError, pd.errors.EmptyDataError):
            report.write("No viral families identified in contigs.\n\n")

print(f"Report markdown successfully generated at {output_md}")
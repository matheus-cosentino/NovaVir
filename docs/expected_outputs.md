# Expected Outputs & Features

## Expected Outputs

Once successfully executed, DeepVir cleanly organizes your results within your standard `--output` directories. Below is a generic overview of an assembled mapping workflow:

```text
results/
├── kraken2_all/                   # Merged Kraken2 outputs and Alpha-Diversity
│   ├── all_samples.biom           # Biom format table encompassing all evaluated samples
│   ├── Rarefaction_Curve.pdf      # Visual plot depicting species richness estimates
│   └── OTU_table.tab              # Tabular format of the operational taxonomic units
├── SRR10677983/
│   ├── log/                       # Sample-specific execution logs
│   ├── fastp/                     # QC & Trimmed adapter reports
│   │   ├── SRR10677983_unp.html                          # Fastp interactive HTML report
│   │   └── SRR10677983_unp.json                          # Fastp metrics in JSON format
│   ├── kraken2_reads/             # Fast Kraken2 taxonomic classification on reads
│   │   ├── SRR10677983_paired_reads_report.txt           # Standard Kraken2 taxonomic report for paired reads
│   │   ├── SRR10677983_paired_reads_output.txt           # Detailed per-read classification output for paired reads
│   │   ├── SRR10677983_paired_reads_biom.txt             # Kraken2 results in BIOM format for paired reads
│   │   ├── SRR10677983_unpaired_reads_report.txt         # Standard Kraken2 taxonomic report for unpaired reads
│   │   ├── SRR10677983_unpaired_reads_output.txt         # Detailed per-read classification output for unpaired reads
│   │   └── SRR10677983_unpaired_reads_biom.txt           # Kraken2 results in BIOM format for unpaired reads
│   ├── diamond_reads/             # Diamond reads classification outputs
│   │   ├── SRR10677983_reads_report.daa                  # DIAMOND alignment in DAA format (intermediate)
│   │   ├── SRR10677983_reads_report.txt                  # Tabular alignment results (outfmt 6)
│   │   └── SRR10677983_reads_hits_with_lineage.tsv       # Diamond hits with full NCBI lineage
│   ├── megan_reads/               # MEGAN Last Common Ancestor taxonomic classification
│   │   └── SRR10677983_reads_summary.megan               # Extracted MEGAN summary table
│   ├── spades/                    # Assembled Fasta contigs
│   │   └── kmer_auto/             # De Novo Assemblage done by SPAdes auto kmer definition           
│   │       └── contigs.fasta                             # Final assembled contigs in FASTA format
│   ├── diamond_spades_kauto/      # Annotated Diamond outputs vs NR
│   │   ├── diamond.log                                   # DIAMOND alignment execution log
│   │   ├── SRR10677983_spades_kauto_report.txt           # Standard alignment report against NR database
│   │   ├── SRR10677983_spades_kauto_hits_with_taxid.tmp  # Intermediate file matching hits with Taxonomic IDs
│   │   ├── SRR10677983_spades_kauto_hits_with_header.tsv # Alignment hits with informative tab-separated headers
│   │   ├── SRR10677983_spades_kauto_hits_no_lineage.temp # Intermediate file before assigning NCBI full lineage
│   │   └── SRR10677983_spades_kauto_hits_with_lineage.tsv # Final comprehensive DIAMOND alignment results with NCBI lineages
│   ├── kraken2_spades_kauto/      # Fast Kraken2 K-mer validations
│   │   ├── SRR10677983_spades_kauto_contig_report.txt    # Standard Kraken2 report for the assembled contigs
│   │   ├── SRR10677983_spades_kauto_contig_output.txt    # Detailed per-contig classification output
│   │   └── SRR10677983_spades_kauto_contig_biom.txt      # Contig taxonomic distribution in BIOM format
│   ├── rvdb_spades_kauto/         # RVDB Summaries containing precise NCBI translation
│   │   ├── SRR10677983_RVDB_results.tsv                  # Primary tabular results evaluated against RVDB
│   │   ├── SRR10677983_RVDB_Summary.csv                  # Aggregated CSV summary report of RVDB profile alignments
│   │   ├── SRR10677983_structural_Orfs.fasta             # ORFs translated from hits to viral structural proteins 
│   │   ├── SRR10677983_structural_contigs.fasta          # Original nucleotide contigs mapped to viral structural proteins
│   │   ├── SRR10677983_pol_Orfs.fasta                    # ORFs translated from hits to viral polymerase proteins
│   │   └── SRR10677983_pol_contigs.fasta                 # Original nucleotide contigs mapped to viral polymerase proteins
│   ├── darkmatter_spades_kauto/   # ORFs & Filtered RdRp viral candidates!
│   │   ├── SRR10677983_ORFs.fasta                        # Multi-table translated ORFs from unmapped contigs
│   │   ├── SRR10677983_RdRp.fev                          # Palmprint raw structural feature evaluation hits
│   │   ├── SRR10677983_RdRp.fasta                        # Sequence data of putative RdRp alignments
│   │   ├── SRR10677983_RdRp.tsv                          # Tabular summary of RdRp predictions
│   │   ├── SRR10677983_Report_Diversity.html             # Comprehensive interactive RMarkdown HTML report
│   │   ├── SRR10677983_contigs_summary.tsv               # Assembly length metrics (Mean, Quartiles, etc.)
│   │   ├── SRR10677983_contigs_summary_per_Kingdom.tsv   # Contig taxonomic breakdown at Kingdom level
│   │   ├── SRR10677983_KINGDOM_Classification_Log10.png  # Barplot comparing Kingdom-level hits
│   │   ├── SRR10677983_Viral_Diamond.tsv                 # Subset of Diamond hits assigned to Viruses
│   │   ├── SRR10677983_Virus_Classification_Log10.png    # Barplot of Viral Family diversity
│   │   ├── SRR10677983_viral_contigs_summary.tsv         # Assembly metrics exclusively for viral families
│   │   ├── SRR10677983_{Family_Name}.fasta               # Disaggregated contigs per viral family (e.g., Flaviviridae)
│   │   ├── SRR10677983_{query}RdRp_Motifs.png            # Palmprint RdRp motif visual diagrams layout per hit
│   │   ├── SRR10677983_RdRp_Orfs.fasta                   # Filtered FASTA containing only target RdRp ORFs
│   │   └── SRR10677983_RdRp_contigs.fasta                # Original nucleotide contigs embedding the RdRp hit
│   └── krona_spades_kauto/        # Interactive HTML plots of relative abundances
│       ├── SRR10677983_spades_kauto_kraken2_krona.html   # Interactive Krona chart for Kraken2 taxonomy
│       └── SRR10677983_spades_kauto_diamond_krona.html   # Interactive Krona chart for DIAMOND taxonomy
├── multiqc_report.html            # Aggregated sequencing QC metrics dashboard
├── multiqc_data/                  # Exported MultiQC data directory
└── stats/                        # Per-sample and per-tool summary statistics for assembly, diamond, kraken2 and dark matter
```




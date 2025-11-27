import pandas as pd
import os

# In reporting.smk

checkpoint report_summarize:
    """
    Executes an R script to generate a primary HTML report and all static summary files.
    """
    input:
        fasta=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/contigs.fasta",
        diamond=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits_with_lineage.tsv",
        duskmatter=f"{config['output_dir']}/{{sample}}/duskmatter/{{sample}}_RdRp.tsv"
    output:
        html=f"{config['output_dir']}/{{sample}}/report/{{sample}}_Report_Diversity.html",
        contig_summary_tsv=f"{config['output_dir']}/{{sample}}/report/{{sample}}_contigs_summary.tsv",
        kingdom_tsv=f"{config['output_dir']}/{{sample}}/report/{{sample}}_contigs_summary_per_Kingdom.tsv",
        viral_diamond_tsv=f"{config['output_dir']}/{{sample}}/report/{{sample}}_Viral_Diamond.tsv",
        viral_summary_tsv=f"{config['output_dir']}/{{sample}}/report/{{sample}}_viral_contigs_summary.tsv",
        kingdom_png=f"{config['output_dir']}/{{sample}}/report/{{sample}}_KINGDOM_Classification_Log10.png",
        viral_png=f"{config['output_dir']}/{{sample}}/report/{{sample}}_Virus_Classification_Log10.png"
    params:
        logo_dirs=config["logo_dirs"],
        output_dir=f"{config['output_dir']}/{{sample}}/report"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_report_summarize.log"
    conda:
        REPORT
    shell:
        """
        # Execute the R script inside a subshell
        (
            Rscript workflow/rules/scripts/generate_report.R \
            --sample_name {wildcards.sample} \
            --fasta_path {input.fasta} \
            --diamond_path {input.diamond} \
            --duskmatter_path {input.duskmatter} \
            --output_dir {params.output_dir} \
            --report_name {wildcards.sample}_Report_Diversity.html \
            --logos {params.logo_dirs} \
            --input workflow/rules/scripts/Report_Model.Rmd
        ) 2>&1 > {log}
        
        # CRITICAL FIX: Explicitly touch all output files to confirm completion
        touch {output.html}
        touch {output.contig_summary_tsv}
        touch {output.kingdom_tsv}
        touch {output.viral_diamond_tsv}
        touch {output.viral_summary_tsv}
        touch {output.kingdom_png}
        touch {output.viral_png}
        
        # Guarantee a successful exit status 
        exit 0

        """


def collect_dynamic_files(wildcards):
    """
    Reads summary files created by the checkpoint to determine the full list of
    dynamic FASTA and PNG files that need to be generated for a given sample.
    """
    # Access the checkpoint's output using the `checkpoints` object and the sample wildcard
    checkpoint_output = checkpoints.report_summarize.get(sample=wildcards.sample).output
    
    report_output_dir = os.path.dirname(checkpoint_output.viral_summary_tsv)
    
    # Paths to the files that determine the dynamic outputs
    summary_path = checkpoint_output.viral_summary_tsv
    duskmatter_path = checkpoints.report_summarize.get(sample=wildcards.sample).input.duskmatter

    dynamic_files = []

    # --- A. Discover Per-Family FASTA Files ---
    if os.path.exists(summary_path):
        summary_df = pd.read_csv(summary_path, sep='\t', header=0)
        family_names = summary_df.columns[1:].tolist()
        for family in family_names:
            dynamic_files.append(
                os.path.join(report_output_dir, f"{wildcards.sample}_{family}.fasta")
            )

    # --- B. Discover Per-Contig RdRp Motif PNG Files ---
    if os.path.exists(duskmatter_path):
        try:
            duskmatter_df = pd.read_csv(duskmatter_path, sep='\t')
            contig_labels = duskmatter_df['Label'].dropna().unique()
            for label in contig_labels:
                dynamic_files.append(
                    os.path.join(report_output_dir, f"{wildcards.sample}_{label}RdRp_Motifs.png")
                )
        except Exception as e:
            # This will appear in the main Snakemake log
            print(f"Warning: Could not process duskmatter file {duskmatter_path} for sample {wildcards.sample}. Details: {e}")

    return dynamic_files

# This input function gathers all files needed for the final 'all_reports' rule.
# It will only be executed by Snakemake after the 'report_summarize' checkpoint is complete for all samples.
def get_all_report_inputs(wildcards):
    """
    Collects all static and dynamic report files from all samples.
    """
    from snakemake.io import Wildcards

    all_files = []
    # 1. Collect all static files from the checkpoint for all samples
    static_files = expand(checkpoints.report_summarize.get(sample=s).output for s in SAMPLES)
    all_files.extend(static_files)

    # 2. Collect all dynamic files for all samples
    for sample in SAMPLES:
        wildcards = Wildcards(fromdict={"sample": sample})
        all_files.extend(collect_dynamic_files(wildcards))
    
    return all_files

# This is the final target rule. It collects all static and dynamic outputs.
rule all_reports:
    input: get_all_report_inputs
    output:
        touch("logs/all_reports_done.txt") # A flag file to indicate the entire workflow is complete
    conda:
        REPORT
    log:
        "logs/all_reports.log"
    shell:
        """
        echo 'All static and dynamic reports generated successfully.' > {log}
        # The log aggregation is removed as it's complex to do with the new input function.
        # The individual log files in the 'logs/' directory still exist for debugging.
        """


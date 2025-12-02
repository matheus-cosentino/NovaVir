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


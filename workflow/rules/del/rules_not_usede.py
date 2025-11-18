# Rules i am not interest anymore
# Cosentino 2025, Novembre 

rule sort_lineage_by_taxonomy:
    """
    Sorts the final output file by the taxonomic lineage columns.
    """
    input:
        "results/diamond/{sample}_{source}_hits_with_lineage.tsv"
    output:
        "results/diamond/{sample}_{source}_hits_with_lineage.sorted.tsv"
    shell:
        """
        (
            # Print the header line first
            head -n 1 {input};
            # Then, sort the rest of the file (from line 2 onwards)
            tail -n +2 {input} | \
            sort -t$'\\t' -k14,14 -k15,15 -k16,16 -k17,17 -k18,18 -k19,19 -k20,20
        ) > {output}
        """

# workflow/rules/summarize_taxonomy.smk

rule summarize_taxonomy_by_rank:
    input:
        hits=f"{config['output_dir']}/diamond/{{sample}}_{{source}}_hits.tsv",
        nodes=config["db"]["taxonnodes"],
        names=config["db"]["taxonnames"],
        accession_map=config["db"]["taxonmap"]
    output:
        summary=f"{config['output_dir']}/summary/{{sample}}_{{source}}_{{rank}}_counts.tsv"
    params:
        rank="{rank}"
    log:
        f"{config['output_dir']}/logs/summary_{{sample}}_{{source}}_{{rank}}.log"
    conda:
        # This environment needs pandas
        SCRIPTS
    script:
        "workflow/rules/scripts/summarize_taxonomy.py"

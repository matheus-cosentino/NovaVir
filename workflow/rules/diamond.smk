# workflow/rule<s/diamond.smk

# This rule for contigs is correct and unchanged.
rule diamond_blastx_contigs:
    input:
        contigs=f"{config['output_dir']}/{{sample}}/assembly/{{sample}}/contigs.fasta"

    output:
        hits=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_contigs_hits.tsv"
    shadow: "minimal" 
    params:
        db=f"{workflow.basedir}/{config['db']['diamond']}",
        #db=config["db"]["diamond"],
        outfmt=config["params"]["diamond"]["outfmt"],
        max_target_seqs=config["params"]["diamond"]["max_target_seqs"],
        evalue=config["params"]["diamond"]["evalue"]
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_contigs.log"
    conda:
        DIAMOND
    shell:
        """
        diamond blastx \
            --query {input.contigs} \
            --db {params.db} \
            --out {output.hits} \
            --threads {resources.threads} \
            --outfmt {params.outfmt} \
            --max-target-seqs {params.max_target_seqs} \
            --evalue {params.evalue} \
            &> {log}
        """

# In workflow/rules/diamond.smk - ROBUST VERSION
rule diamond_blastx_reads:
    input:
        r1=f"{config['output_dir']}/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{config['output_dir']}/trimmed/{{sample}}_2.fastq.gz",
        orphans=f"{config['output_dir']}/trimmed/{{sample}}_orphans.fastq.gz",
        layout=f"{config['output_dir']}/sra_layout/{{sample}}.layout.txt"
    output:
        hits=f"{config['output_dir']}/{{sample}}/diamond/{{sample}}_reads_hits.tsv"
    shadow: 
        "minimal" 
    params:
        db=f"{workflow.basedir}/{config['db']['diamond']}",
        outfmt=config["params"]["diamond"]["outfmt"],
        max_target_seqs=config["params"]["diamond"]["max_target_seqs"],
        evalue=config["params"]["diamond"]["evalue"]
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_reads.log"
    conda:
        DIAMOND
    shell:
        """
        # Read the layout from the input file
        LAYOUT=$(cat {input.layout})

        # Build the query file string based on the layout
        if [ "$LAYOUT" = "PE" ]; then
            QUERY_FILES="--query {input.r1} --query {input.r2} --query {input.orphans}"
        else
            QUERY_FILES="--query {input.r1}"
        fi

        diamond blastx \
            $QUERY_FILES \
            --db {params.db} \
            --out {output.hits} \
            --threads {resources.threads} \
            --outfmt {config[params][diamond][outfmt]} \
            --max-target-seqs {config[params][diamond][max_target_seqs]} \
            --evalue {config[params][diamond][evalue]} \
            &> {log}
        """

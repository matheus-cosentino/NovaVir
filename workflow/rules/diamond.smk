###################################################################################
#                       workflow/rules/diamond.smk                                #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                                                                                 #
# oooooooooo.    o8o                               oooooo     oooo  o8o           #
# `888'   `Y8b   `"'                                `888.     .8'   `"'           #
#  888      888 oooo   .oooo.o  .ooooo.   .ooooo.    `888.   .8'   oooo  oooo d8b #
#  888      888 `888  d88(  "8 d88' `"Y8 d88' `88b    `888. .8'    `888  `888""8P #
#  888      888  888  `"Y88b.  888       888   888     `888.8'      888   888     #
#  888     d88'  888  o.  )88b 888   .o8 888   888      `888'       888   888     #
# o888bood8P'   o888o 8""888P' `Y8bod8P' `Y8bod8P'       `8'       o888o d888b    #
#                                                                                 #
###################################################################################
#                              version: 12.2025                                   #
###################################################################################

rule diamond_blastx_contigs:
    input:
       contigs = get_contigs_path
    output:
        hits="{out_dir}/{sample}/diamond_{tool}/{sample}_contigs_report.txt"
    params:
        #db=f"{workflow.basedir}/{config['resources']['diamond']}",
        db=config["resources"]["diamond"],
        outfmt=config["diamond"]["outfmt"],
        max_target_seqs=config["diamond"]["max_target_seqs"],
        evalue=config["diamond"]["evalue"]
    log:
        "{out_dir}/{sample}/log/diamond_contigs_{tool}_{sample}.log"
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
            --log \
            &> {log}
        """








#develop
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
            --log \
            &> {log}
        """

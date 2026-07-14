###################################################################################
#                      workflow/rules/reportin.smk                                #
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
#                              version: 03.2026                                   #
###################################################################################

#######################################################################
# --- Statistics Generation for MultiQC Integration ---              #
#######################################################################

# --- Generate Assembly Statistics ---
rule generate_assembly_stats:
    """Generate contig statistics from assembly output"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        fasta = get_contigs_path
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_assembly_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_assembly_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_assembly_stats_{tool}_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_assembly_stats.py"


# --- Generate Diamond Statistics ---
rule generate_diamond_stats:
    """Generate Diamond hit statistics for contigs"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        diamond = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_diamond_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_diamond_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_diamond_stats_{tool}_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_diamond_stats.py"


# --- Generate Diamond Statistics (Reads) ---
rule generate_diamond_reads_stats:
    """Generate Diamond hit statistics for reads"""
    input:
        diamond = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_reads_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_reads_diamond_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_reads_diamond_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_diamond_reads_stats_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_diamond_stats.py"


# --- Generate Kraken2 Statistics (Contigs) ---
rule generate_kraken2_stats:
    """Generate Kraken2 classification statistics for contigs"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_kraken2_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_kraken2_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_kraken2_stats_{tool}_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_kraken2_stats.py"


# --- Generate Kraken2 Statistics (Reads) ---
rule generate_kraken2_reads_stats:
    """Generate Kraken2 classification statistics for reads"""
    input:
        report = lambda wc: os.path.join(OUT_DIR, wc.sample, "kraken2_reads", 
                            f"{wc.sample}_{'paired' if (SAMPLE_META.get(wc.sample, {}).get('mode') == 'PAIRED' or (SAMPLE_META.get(wc.sample, {}).get('mode') == 'SRA' and len(SAMPLE_META.get(wc.sample, {}).get('files', [])) == 2)) else 'unpaired'}_reads_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_reads_kraken2_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_reads_kraken2_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_kraken2_reads_stats_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_kraken2_stats.py"


# --- Generate Dark Matter Statistics ---
rule generate_darkmatter_stats:
    """Generate dark matter (no-hit contigs) statistics"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        nohits = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_{tool}_nohit.fasta")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_darkmatter_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_{tool}_darkmatter_stats.tsv")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "generate_darkmatter_stats_{tool}_{sample}.log")
    conda:
        CORE
    script:
        "../scripts/generate_darkmatter_stats.py"

rule multiqc_aggregate:
    message:
        """
        Generate a multiqc report in HTML format all samples
        """
    conda:
        MULTIQC
    input:
        files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)],
        final_outputs = get_final_outputs()
    output:
        report = os.path.join(OUT_DIR, "multiqc_report.html"),
        data_dir = directory(os.path.join(OUT_DIR, "multiqc_data"))
    params:
        files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)],
        #config_override = "sp: { diamond/log: { fn: '*_diamond.log' } }",
        extra = "--title 'DiscoVir Aggregate Report'"
    conda:
        MULTIQC    
    log:
        os.path.join(OUT_DIR, "log", "multiqc_agregate.log")
    shell:
        """
        multiqc \
        --quiet \
        --export \
        --force \
        --outdir {OUT_DIR} \
        --filename multiqc_report.html \
        {params.extra} \
        {input.files} > {log} 2>&1 
        """

# --- Rule 2: Per-Sample Report ---
rule multiqc_sample:
    message:
        "Generate MultiQC report for sample: {wildcards.sample}"
    conda:
        MULTIQC
    input:
        # We pass the wildcards object to get inputs ONLY for this sample
        files = lambda wc: get_multiqc_inputs(wildcards=wc)
    output:
        report = os.path.join(OUT_DIR, "{sample}", "multiqc", "multiqc_report.html"),
        data_dir = directory(os.path.join(OUT_DIR, "{sample}", "multiqc"))
    params:
        #config_override = "sp: {{ diamond/log: {{ fn: '*_diamond.log' }} }}",
        extra = "--title 'Report for {sample}'"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "multiqc.log")
    shell:
        """
        multiqc \
        --quiet \
        --export \
        --force \
        --outdir {output.data_dir} \
        --filename multiqc_report.html \
        {params.extra} \
        {input.files} > {log} 2>&1 
        """
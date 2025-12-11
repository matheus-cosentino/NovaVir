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
#                              version: 12.2025                                   #
###################################################################################
#
#rule snakemake_report:
    # Aim: generates a workflow report in HTML format
    # Use: snakemake --report [OPTIONS] [REPORT]
#    message:
#        """
#        Generate a workflow report in HTML format 
#        """
#    conda:
#        CORE
#    input:
#        final_outputs = get_final_outputs()
#    output:
#        html_report = os.path.join(OUT_DIR, "{sample}", "reporting", "{sample}_reporting.html")
#    log:
#        os.path.join(OUT_DIR, "{sample}", "log", "reporting_{sample}.log")
#    shell:
#        "snakemake "            # Snakemake
#        "--report "              # Create an HTML report with results and statistics
#        " {output.html_report} " # Output report
#        "2> {log}"               # Log redirection


rule multiqc_aggregate:
    message:
        """
        Generate a multiqc report in HTML format all samples
        """
    conda:
        MULTIQC
    input:
        final_outputs = get_final_outputs()
        #files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)]
    output:
        report = os.path.join(OUT_DIR, "multiqc_all", "multiqc_report.html"),
        data_dir = directory(os.path.join(OUT_DIR, "multiqc_all"))
    params:
        files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)],
        config_override = "sp: { diamond/log: { fn: '*_diamond.log' } }",
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
        --outdir {output.data_dir} \
        --filename {output.report} \
        --cl-config "{params.config_override}" \
        {params.extra} \
        {params.files} > {log} 2>&1 
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
        config_override = "sp: {{ diamond/log: {{ fn: '*_diamond.log' }} }}",
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
        --filename {output.report} \
        --cl-config "{params.config_override}" \
        {params.extra} \
        {input.files} > {log} 2>&1 
        """
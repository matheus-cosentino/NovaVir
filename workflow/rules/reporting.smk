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


rule multiqc_report:
    message:
        """
        Generate a multiqc report in HTML format 
        """
    conda:
        MULTIQC
    input:
        final_outputs = get_final_outputs() ,
        files = get_multiqc_inputs()
    output:
        report = os.path.join(OUT_DIR, "multiqc", "multiqc_report.html"),
        data_dir = directory(os.path.join(OUT_DIR, "multiqc"))
    params:
        config_override = "sp: { diamond/log: { fn: '*_diamond.log' } }"
    conda:
        MULTIQC    
    log:
        os.path.join(OUT_DIR, "log", "multiqc.log")
    shell:
        """
        multiqc \
        --quiet \
        --export \
        --outdir {output.data_dir} \
        --cl-config "{params.config_override}" \
        {input.files} > {log} 2>&1 """
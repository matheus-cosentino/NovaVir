###################################################################################
#                        workflow/rules/fastp.smk                                 #
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



############################################
# --- 1. Filter fastq.gz Paired Files --- #
########################################### 

#rule to process paired data
rule fastp_paired:
    input:
        r1 = get_input_r1,
        r2 = get_input_r2
    output:
        # Define all possible output files explicitly and statically.
        #r1= temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_1.fastq.gz")),
        #r2= temp(os.path.join(OUT_DIR, "{sample}", "fastp" , "{sample}_2.fastq.gz")),
        r1= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_1.fastq.gz"),
        r2= os.path.join(OUT_DIR, "{sample}", "fastp" , "{sample}_2.fastq.gz"),
        orphans=temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_orphans.fastq.gz")),
        html = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_paired.html"),
        json= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_paired.json")
    log:
        #f"{config['output_dir']}/{{sample}}/logs/{{sample}}_fastp.log"
        #"{out_dir}/{sample}/log/{sample}_fastp_paired.log"
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_fastp_paired.log")
    params:
        length_required=config['fastp']['length_required'],
        quality=config['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --unpaired1 {output.orphans} --unpaired2 {output.orphans} \
            --html {output.html} \
            --json {output.json} \
            --thread {resources.threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} \
            2> {log}
        """

#if fail to make the paired reads cleanup, follown immediatelly to the fastp_unpaired rule
#############################################
# --- 2. Filter fastq.gz unpaired Files --- #
#############################################

rule fastp_unpaired:
    input:
        reads = get_input_unp,
    output:
        #r1= temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.fastq.gz")),
        r1= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.fastq.gz"),
        html= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.html"),      
        json= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.json")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_fastp_unpaired.log")
    params:
        length_required=config['fastp']['length_required'],
        quality=config['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.reads} --out1 {output.r1} \
            --html {output.html} \
            --json {output.json} \
            --thread {resources.threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} \
            --dedup --trim_poly_g --trim_poly_x \
            2> {log}
        """

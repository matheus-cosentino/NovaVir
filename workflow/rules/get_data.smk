# workflow/rules/get_data.smk

checkpoint determine_sra_layout:
    """
    Checks downloaded files to determine if a sample is PE or SE.
    It works by checking if the R2 file has a size greater than zero.
    """
    input:
        # This checkpoint now depends on the output of download_sra_data
        r1="raw_data/{sample}_1.fastq.gz",
        r2="raw_data/{sample}_2.fastq.gz"
    output:
        # The output is the same: a file containing either "PE" or "SE"
        f"{config['output_dir']}/{{sample}}/sra_layout/{{sample}}.layout.txt"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_sra_layout.log"
    shell:
        """
        # The '-s' test in shell checks if a file exists AND has a size greater than zero.
        # This is the perfect way to distinguish a real R2 file from an empty placeholder.
        if [ -s {input.r2} ]; then
            echo "PE" > {output}
        else
            echo "SE" > {output}
        fi
        """


rule download_sra_data:
    output:
        r1="raw_data/{sample}_1.fastq.gz",
        r2="raw_data/{sample}_2.fastq.gz"
    log:
        f"{config['output_dir']}/{{sample}}/logs/{{sample}}_download.log"
    conda:
        DOWNLOAD
    params:
        # Wrap the string in a lambda function to delay wildcard evaluation.
        tmpdir=lambda wildcards: f"raw_data/tmp_{wildcards.sample}"
    shadow: 
        "minimal" 
    shell:
        """
        # Create the temporary directory defined in params
        mkdir -p {params.tmpdir}

        # Run fasterq-dump, outputting to the temporary directory
        fasterq-dump --threads {resources.threads} --split-files -O {params.tmpdir} {wildcards.sample} > {log} 2>&1

        # Gzip the R1 file and move it
        gzip {params.tmpdir}/{wildcards.sample}_1.fastq
        mv {params.tmpdir}/{wildcards.sample}_1.fastq.gz {output.r1}

        # Check if the R2 file was created
        if [ -f {params.tmpdir}/{wildcards.sample}_2.fastq ]; then
            # Paired-End: Gzip and move the R2 file.
            gzip {params.tmpdir}/{wildcards.sample}_2.fastq
            mv {params.tmpdir}/{wildcards.sample}_2.fastq.gz {output.r2}
        else
            # Single-End: Create an empty placeholder R2 file.
            touch {output.r2}
        fi
        
        # Clean up the temporary directory
        rm -r {params.tmpdir}
        """

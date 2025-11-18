# workflow/rules/reporting.smk

rule generate_markdown_report:
    input:
        fastp_reports=expand(f"{config['output_dir']}/reports/fastp/{{sample}}.json", sample=SAMPLES),
        spades_contigs=expand(f"{config['output_dir']}/assembly/{{sample}}/contigs.fasta", sample=SAMPLES),
        read_hits=expand(f"{config['output_dir']}/diamond/{{sample}}_reads_hits.tsv", sample=SAMPLES),
        contig_hits=expand(f"{config['output_dir']}/diamond/{{sample}}_contigs_hits.tsv", sample=SAMPLES),
        read_summaries=expand(f"{config['output_dir']}/summary/{{sample}}_reads_{config['params']['summary']['rank']}_counts.tsv", sample=SAMPLES),
        contig_summaries=expand(f"{config['output_dir']}/summary/{{sample}}_contigs_{config['params']['summary']['rank']}_counts.tsv", sample=SAMPLES),
        logo=config["report"]["logo"]
    output:
        report_md=f"{config['output_dir']}/final_summary_report.md"
    params:
        title=config["report"]["title"],
        author=config["report"]["author"]
    log:
        "logs/final_report.log"
    conda:
        SCRIPTS
    script:
        workflow.source_path("scripts/generate_report.py")

rule compile_pdf_report:
    input:
        report_md=f"{config['output_dir']}/final_summary_report.md",
        logo=config["report"]["logo"]
    output:
        report_pdf=f"{config['output_dir']}/final_summary_report.pdf"
    log:
        "logs/pandoc_report.log"
    conda:
        PANDOC
    shell:
        """
        pandoc {input.report_md} \
            -o {output.report_pdf} \
            --pdf-engine=xelatex \
            -V geometry:"margin=1in" \
            &> {log}
        """
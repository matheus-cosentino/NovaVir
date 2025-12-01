#!/usr/bin/env Rscript

# Este script é um wrapper que usa a diretiva 'script:' do Snakemake
# para chamar a função GenerateSinglePlot do script 'generate_dynamic_outputs.R'.

# --- 1. Capturar Parâmetros do Snakemake ---
duskmatter_path <- snakemake@input$duskmatter
contig_label    <- snakemake@wildcards$label
sample_name     <- snakemake@wildcards$sample
output_path     <- snakemake@output[[1]]
log_file        <- snakemake@log[[1]]

# O caminho para o script que contém a função
source_script   <- snakemake@params$dynamic_outputs_script

# --- 2. Executar a Função e Capturar a Saída ---

# Redireciona toda a saída (stdout e stderr) para o arquivo de log
sink(log_file, append = TRUE, type = "output")
sink(log_file, append = TRUE, type = "message")

source(source_script) # Carrega o script com a função GenerateSinglePlot
GenerateSinglePlot(duskmatter_path, contig_label, sample_name, output_path)

sink(type = "message"); sink() # Restaura a saída padrão
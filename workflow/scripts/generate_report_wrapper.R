#!/usr/bin/env Rscript

# Este script é um wrapper que usa a diretiva 'script:' do Snakemake
# para chamar o script 'generate_report.R', que espera argumentos de linha de comando.

# --- 1. Verificação de Pacotes ---
if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
library(stringr)

# --- 2. Capturar Parâmetros do Snakemake ---
# O Snakemake disponibiliza o objeto 'snakemake' automaticamente
sample_name       <- snakemake@wildcards$sample
fasta_path        <- snakemake@input$fasta
diamond_path      <- snakemake@input$diamond
duskmatter_path   <- snakemake@input$duskmatter
output_dir        <- snakemake@params$output_dir
report_template   <- snakemake@params$report_template
logo_dirs         <- snakemake@params$logo_dirs
log_file          <- snakemake@log[[1]]

# O nome do relatório é derivado da saída 'html' e o script a ser chamado é pego dos params
report_name       <- basename(snakemake@output$html)
script_a_chamar   <- snakemake@params$report_script

# --- 3. Construir e Executar o Comando ---

# Construir a string de argumentos
args <- c(
    "--sample_name", sample_name,
    "--fasta_path", fasta_path,
    "--diamond_path", diamond_path,
    "--duskmatter_path", duskmatter_path,
    "--output_dir", output_dir,
    "--report_name", report_name,
    "--logos", paste0('"', logo_dirs, '"'), # Adiciona aspas para segurança
    "--input", report_template
)

# Executar o script original, redirecionando a saída para o arquivo de log
system2("Rscript", args = c(script_a_chamar, args), stdout = log_file, stderr = log_file)

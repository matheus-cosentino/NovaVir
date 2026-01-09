#!/bin/bash
#SBATCH --job-name=DiscoVir
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs_slurm/DiscoVir_%j.out

echo "
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
##                                                      ##
##      ▗▄▄▄ ▗▄▄▄▖ ▗▄▄▖ ▗▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄        ##
##      ▐▌  █  █  ▐▌   ▐▌   ▐▌ ▐▌▐▌  ▐▌  █  ▐▌ ▐▌       ##
##      ▐▌  █  █   ▝▀▚▖▐▌   ▐▌ ▐▌▐▌  ▐▌  █  ▐▛▀▚▖       ##
##      ▐▙▄▄▀▗▄█▄▖▗▄▄▞▘▝▚▄▄▖▝▚▄▞▘ ▝▚▞▘ ▗▄█▄▖▐▌ ▐▌       ##
##                                                      ##
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
## Name ___________________ DiscoVir.slurm              ##
## Version ________________ v.2025.12                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ PGGEN ##
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
"

# --- Definições de Diretórios ---
SUBMIT_DIR=$SLURM_SUBMIT_DIR
SCRATCH_DIR="/scratch/$USER/DiscoVir_${SLURM_JOB_ID}"

echo "----------------------------------------------------------------"
echo "Job ID: $SLURM_JOB_ID"
echo "Submit Dir: $SUBMIT_DIR"
echo "Scratch Dir: $SCRATCH_DIR"
echo "----------------------------------------------------------------"

# --- 1. Setting Scratch ---
mkdir -p "$SCRATCH_DIR"
mkdir -p "$SCRATCH_DIR/logs_slurm"

echo "[INFO] Sending Files to Scratch with rsync..."

rsync -av --exclude '.snakemake' \
    "$SUBMIT_DIR/config" \
    "$SUBMIT_DIR/workflow" \
    "$SUBMIT_DIR/resources" \
    "$SUBMIT_DIR/DiscoVir.sh" \
    "$SUBMIT_DIR/profile_slurm" \
    "$SUBMIT_DIR/accessions.txt" \
    "$SCRATCH_DIR/"

# Link simbólico para os dados (evita copiar gigabytes de dados brutos)
# O script assume que seus dados estão em 'data' no diretório de submissão
ln -sfn "$SUBMIT_DIR/data" "$SCRATCH_DIR/data"

# --- 2. Configurar Ambiente ---
source /usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh

# --- 3. Executar o Pipeline no Scratch ---
cd "$SCRATCH_DIR"

echo "[INFO] Starting DiscoVir..."

# Running DiscoVir. 
bash DiscoVir.sh \
    --input data \
    --output results \
    --sra accessions.txt \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db /share/banks/diamond/nr/21-09/diamond_nr.dmnd \
    --kraken2_db /share/banks/kraken2/nt/21-09/nt/

EXIT_CODE=$?

# --- 4. Retrieving Data ---
echo "[INFO] Trazendo resultados de volta para o Home com rsync..."

# Define onde os resultados finais ficarão no seu Home
FINAL_OUTPUT_DIR="$SUBMIT_DIR/results_job_${SLURM_JOB_ID}"
mkdir -p "$FINAL_OUTPUT_DIR"

# O rsync aqui é crucial.
# -a: archive mode (preserva permissões, donos, grupos, tempos)
# -v: verbose
# Copiamos apenas a pasta 'results'. 
# NÃO trazemos a pasta '.snakemake' (que contém os temp files pesados do spades que causaram o erro 28),
# a menos que você precise debugar logs específicos.
rsync -av "$SCRATCH_DIR/results/" "$FINAL_OUTPUT_DIR/"

# Opcional: Trazer logs do snakemake para debug se falhar
if [ $EXIT_CODE -ne 0 ]; then
    echo "[WARN] O Job falhou. Trazendo logs do Snakemake para análise."
    mkdir -p "$FINAL_OUTPUT_DIR/snakemake_logs"
    rsync -av "$SCRATCH_DIR/.snakemake/log/" "$FINAL_OUTPUT_DIR/snakemake_logs/"
fi

# Limpeza do Scratch
# echo "[INFO] Limpando Scratch..."
# rm -rf "$SCRATCH_DIR"

echo "----------------------------------------------------------------"
echo "Job finalizado com código: $EXIT_CODE"
echo "Resultados em: $FINAL_OUTPUT_DIR"
exit $EXIT_CODE
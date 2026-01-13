#!/bin/bash
#SBATCH --job-name=DiscoVir_Orch
#SBATCH --output=logs_slurm/orchestrator_%j.log
#SBATCH --error=logs_slurm/orchestrator_%j.err
#SBATCH --time=7-00:00:00        
#SBATCH --mem=2G                 
#SBATCH --cpus-per-task=1
#SBATCH --partition=long         

# --- 1. ENVIRONMENT SETUP ---
echo "[INFO] Starting Orchestrator on host: $(hostname)"
source /usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate DiscoVir

# --- 2. SAFE DIRECTORY DEFINITIONS ---
PROJECT_HOME="$SLURM_SUBMIT_DIR"

# Output Final (Continua no HOME/Storage compartilhado obrigatóriamente)
RESULTS_DIR="$PROJECT_HOME/results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

# [NOVO] Definição do Temp Local do Node
# Se o Slurm definir um TMPDIR, usa ele. Senão, usa /tmp.
NODE_TMP="${SLURM_TMPDIR:-/tmp}"
echo "[INFO] Temp Dir for Nodes set to: $NODE_TMP"

# --- 3. PIPELINE EXECUTION ---
cd "$PROJECT_HOME" || { echo "Error: Failed to access $PROJECT_HOME"; exit 1; }

# Adicionamos a flag --temp-dir
bash DiscoVir.sh \
    --input /home/cosentino/Hanta \
    --output "$RESULTS_DIR" \
    --temp-dir "$NODE_TMP" \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db /databases/tools/diamond/ncbi-nr.taxonomy.2.1.13.dmnd \
    --kraken2_db /databases/tools/kraken2/k2_standard_20250714

if [ $? -eq 0 ]; then
    echo "✅ Pipeline Finished Successfully. Results in: $RESULTS_DIR"
else
    echo "❌ Pipeline failed."
    exit 1
fi
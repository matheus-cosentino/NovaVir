#!/bin/bash
#SBATCH --job-name=DiscoVir
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs_slurm/DiscoVir_%j.out

echo "
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
## Name ___________________ DiscoVir.slurm              ##
## Version ________________ v.2026.07                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ LDDV  ##
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
"

echo "##### DiscoVir SLURM Wrapper #####"

# -------- KEEP ---------
# --- 1. DEFINE PATHS ---
# PROJECT_HOME: Where your Snakefile and Config live (MUST be shared/home)
PROJECT_HOME=${SLURM_SUBMIT_DIR}

#Dynamic Variables per run
DATA_DIR=$1
OUT_DIR=$2

# Fixed Variables
DIAMOND_BANK=/databases/tools/diamond/ncbi-nr.taxonomy.2.1.13.dmnd
KRAKEN2_BANK=/databases/tools/kraken2/k2_standard_20250714

# SCRATCH_DIR: Where the final results go. 
SCRATCH_DIR="/scratch/USERS/$USER/DiscoVir_$SLURM_JOB_ID"
mkdir -p "$SCRATCH_DIR"

# TEMP_DIR: Where TEMP files are generated
TEMP_DIR="/scratch/USERS/$USER/temp_DiscoVir_$SLURM_JOB_ID"

# --- 2. ACTIVATE ENVIRONMENT ---
### load the necessary modules
module load miniconda3
### load the configurations for the conda environment
eval "$(conda shell.bash hook)"

cd "$PROJECT_HOME" || { echo "Error: Could not find $PROJECT_HOME"; exit 1; }
echo "[INFO] Running from Shared Path: $(pwd)"
echo "[INFO] Heavy Data will be written to: $SCRATCH_DIR"

# -------- CHANGE ---------
# --- 3. EXECUTE FROM HOME ---
bash DiscoVir.sh \
    --input ${DATA_DIR} \
    --output ${SCRATCH_DIR} \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db ${DIAMOND_BANK} \
    --kraken2_db ${KRAKEN2_BANK}\
    --tempd-dir ${TEMP_DIR}

# --- 5. RETRIEVE RESULTS ---
# Move results from Scratch back to Home (or a safe place) before the job ends
FINAL_HOME_DIR="/home/${USER}/results/DiscoVir/results_${SLURM_JOB_ID}"
mkdir -p "${FINAL_HOME_DIR}"

echo "[INFO] Retrieving results from Scratch..."
rsync -av "${SCRATCH_DIR}/" "${FINAL_HOME_DIR}/"

echo "[INFO] Done."
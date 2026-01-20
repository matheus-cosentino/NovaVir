#!/bin/bash
#SBATCH --job-name=DiscoVir
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs_slurm/DiscoVir_%j.out
#SBATCH --constraint infiniband

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

echo "##### DiscoVir SLURM Wrapper - Split Storage Mode #####"

# --- 1. DEFINE PATHS ---
# PROJECT_HOME: Where your Snakefile and Config live (MUST be shared/home)
PROJECT_HOME="$SLURM_SUBMIT_DIR"

# SCRATCH_OUT: Where the final results go. 
SCRATCH_DIR="/scratch-ib/$USER/DiscoVir_Results"
mkdir -p "$SCRATCH_DIR"

# --- 2. ACTIVATE ENVIRONMENT ---
source /usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh
conda activate DiscoVir

# --- 3. EXECUTE FROM HOME (The Fix) ---
cd "$PROJECT_HOME" || { echo "Error: Could not find $PROJECT_HOME"; exit 1; }

echo "[INFO] Running from Shared Path: $(pwd)"
echo "[INFO] Heavy Data will be written to: $SCRATCH_DIR"

# --- 4. RUN PIPELINE ---
# We pass the SCRATCH_DIR to the script
bash DiscoVir.sh \
    --input /home/cosentino/discovir/data \
    --output "$SCRATCH_DIR" \
    --assembly \
#    --diamond \
#    --reads-diamond \
#    --reads-kraken \
#    --kraken2 \
    --darkmatter \
#   --diamond_db /databases/tools/diamond/ncbi-nr.taxonomy.2.1.13.dmnd
#   --kraken2_db /databases/tools/kraken2/k2_standard_20250714

if [ $? -eq 0 ]; then
    echo "✅ Pipeline Finished Successfully. Results in: $RESULTS_DIR"
else
    echo "❌ Pipeline failed."
    exit 1
fi
# --- 5. RETRIEVE RESULTS ---
# Move results from Scratch back to Home (or a safe place) before the job ends
FINAL_HOME_DIR="$PROJECT_HOME/results_$SLURM_JOB_ID"
mkdir -p "$FINAL_HOME_DIR"

echo "[INFO] Retrieving results from Scratch..."
rsync -av "$SCRATCH_DIR/" "$FINAL_HOME_DIR/"

echo "[INFO] Done."
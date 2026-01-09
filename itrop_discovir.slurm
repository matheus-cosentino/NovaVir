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
# WARNING: If scratch is local, you might lose results unless this path is shared.
# If your /home is full, try to output to a shared /data or /projects folder if available.
# For now, we will output to /scratch but you must check if results persist.
SCRATCH_DIR="/scratch-ib/$USER/DiscoVir_Results_$SLURM_JOB_ID"
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
    --input data \
    --output "$SCRATCH_DIR" \
    --sra accessions.txt \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db /share/banks/diamond/nr/21-09/diamond_nr.dmnd \
    --kraken2_db /share/banks/kraken2/nt/21-09/nt/

# --- 5. RETRIEVE RESULTS ---
# Move results from Scratch back to Home (or a safe place) before the job ends
FINAL_HOME_DIR="$PROJECT_HOME/results_$SLURM_JOB_ID"
mkdir -p "$FINAL_HOME_DIR"

echo "[INFO] Retrieving results from Scratch..."
rsync -av "$SCRATCH_DIR/" "$FINAL_HOME_DIR/"

echo "[INFO] Done."
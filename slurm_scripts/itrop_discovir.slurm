#!/bin/bash
#SBATCH --job-name=DiscoVir
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs_slurm/DiscoVir_%j.out
#SBATCH --constraint=infiniband

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
## Version ________________ v.2026.03                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ PGGEN ##
#####M#A#C#C##2#0#2#5###############I#R#D##&##U#F#R#J#####
"

echo "##### DiscoVir SLURM Wrapper - Itrop Model #####"

bash DiscoVir.sh \
    --input /home/cosentino/phd/Carollia \
    --output /scratch-ib/$USER/DiscoVir_Demo \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db /databases/tools/diamond/ncbi-nr.taxonomy.2.1.13.dmnd \
    --kraken2_db /databases/tools/kraken2/k2_standard_20250714

echo "[INFO] Done."
#!/bin/bash
###################################################################################
#                               DiscoVir.sh                                       #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                                                                                 #
# oooooooooo.    o8o                               oooooo     oooo  o8o           #
# `888'   `Y8b   `"'                                `888.     .8'   `"'           #
#  888      888 oooo   .oooo.o  .ooooo.   .ooooo.    `888.   .8'   oooo  oooo d8b #
#  888      888 `888  d88(  "8 d88' `"Y8 d88' `88b    `888. .8'    `888  `888""8P #
#  888      888  888  `"Y88b.  888       888   888     `888.8'      888   888     #
#  888     d88'  888  o.  )88b 888   .o8 888   888      `888'       888   888     #
# o888bood8P'   o888o 8""888P' `Y8bod8P' `Y8bod8P'       `8'       o888o d888b    #
#                                                                                 #
###################################################################################
#                              version: 03.2026                                   #
###################################################################################

# --- Color Palettes ---
blu="\033[1;34m"  
green="\033[1;32m" 
red="\033[1;31m"   
ylo="\033[1;33m"   
nc="\033[0m"       

# Nome do ambiente Conda (deve bater com workflow/envs/DiscoVir.yaml)
ENV_NAME="DiscoVir"
ENV_FILE="workflow/envs/DiscoVir.yaml"

#####################
# --- Functions --- #
#####################

# --- Spinner (Safe Version) --- #
run_with_spinner() {
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null
    local sequence="GGACTATCCTAAGTGTGACCGTGCTTTGCCGAGCATGATTAGGATGATTTCTGCCATGATACTTGGCTCTAAGCACACAACTTGCTGCACAAATAGTGATAGGTATTACAGATTGTGCAATGAGTTGGCACAAGTGCTCACTGAAGTTGTTTATTCCAATGGTGGTTTTTATTTTAAACCAGGAGGTACAACTTCAGGTGATGCAACTACAGCATATGCCAATTCTGTTTTCAACATATTCCAGGCTGTCAGTGCTAACATTAACCGTTTGCTCACTGTTGACAGTTATGCTATTCATAATGATTCTGTCAAGAGTTTGCAGAGGCAGTTGTATGACAATTGCTACCGTGCCACTTCTGTA"
    local seq_len=${#sequence}
    local width=15
    local pos=0
    
    # REMOVED: tput civis (caused the crash)

    local cA="\033[1;32m"
    local cT="\033[1;31m"
    local cC="\033[1;34m"
    local cG="\033[1;33m"
    local nc="\033[0m"

    while kill -0 $pid 2>/dev/null; do
        local chunk=""
        for (( j=0; j<width; j++ )); do
            local idx=$(( (pos + j) % seq_len ))
            local base="${sequence:$idx:1}"
            case "$base" in
                A) chunk="${chunk}${cA}A${nc}" ;;
                T) chunk="${chunk}${cT}T${nc}" ;;
                C) chunk="${chunk}${cC}C${nc}" ;;
                G) chunk="${chunk}${cG}G${nc}" ;;
                *) chunk="${chunk}${base}" ;;
            esac
        done
        printf "\r\033[K [ ${chunk} ] Processing..."
        pos=$(( (pos + 1) % seq_len ))
        sleep 0.1
    done
    wait $pid
    local exit_code=$?
    
    # REMOVED: tput cnorm (caused the crash)
    
    printf "\r\033[K"
    if [ $exit_code -eq 0 ]; then
        echo -e "${green}✔ Job done!${nc}"
    else
        echo -e "${red}✖ Job failed with exit code $exit_code.${nc}"
        exit $exit_code
    fi
}

# --- Help --- #
help(){
 echo -e "
 ${green}
 DiscoVir${nc}: Viral Metagenomics & 'Dark Matter' Discovery

 ${green}Author${nc}: MSc. Matheus Cosentino 
 ${green}Version${nc}: 03.2026

 ${ylo}Usage:${nc}
  bash DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 ${ylo}Required Arguments:${nc}
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved
 
 ${ylo}Module Toggles (Enable/Disable Analysis):${nc}
   --assembly           Enable De Novo Assembly (Default: Disabled)
   --kraken2            Enable Kraken2 Contigs Taxonomy (Default: Disabled)
   --diamond            Enable Diamond ContigsTaxonomy (Default: Disabled)
   --darkmatter         Enable Palm Annot / Dark Matter (Default: Disabled)
   --rvdb               Enable RVDB validation on dark matter ORFs (Default: Disabled)
   --download-only      Only download data from SRA provided in --sra list (Default: Disabled)
   --remove-download    Remove downloaded data from SRA (Default: Disabled)
   --reads-kraken       Enable Kraken2 analysis on reads (Default: Disabled)
   --reads-diamond      Enable Diamond analysis on reads, including MEGAN LCA (Default: Disabled)
 
 ${ylo}Database Overrides (Define External Databases):${nc}
   --diamond_db <FILE>  Diamond database file (.dmnd).
   --kraken2_db <DIR>   Kraken2 database directory (Must contain hash.k2d, opts.k2d, taxo.k2d)

 ${ylo}Optional Arguments:${nc}
   -h, --help           Show this help message
   -v, --version        Show version
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of total submitted jobs at same time (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)
   --temp-dir <DIR>     Scratch directory for temporary files AND Snakemake workdir
                        (default: auto-discovered from profile's default-resources.tmpdir)
   --workdir <DIR>      Override only the Snakemake working directory if you need it
                        separate from --temp-dir (default: same as --temp-dir)
 "
}

version(){
    echo "DiscoVir v.03.2026"
}

###################################
# --- Environment Management --- #
###################################

manage_environment(){
    echo -e "${blu}[INFO]${nc} Checking Conda environment '${ylo}${ENV_NAME}${nc}'..."
    
    # 1. Try to run conda hook directly
    if __conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"; then
        eval "$__conda_setup"
    # 2. If that fails (else), use the fallback paths
    else
        if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/miniconda3/etc/profile.d/conda.sh"
        elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/anaconda3/etc/profile.d/conda.sh"
        elif [ -f "/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh" ]; then
            . "/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh"
        else
            echo -e "${red}[ERROR]${nc} Could not find conda.sh. Ensure Conda is installed."
            exit 1
        fi
    fi

    if [[ ! -f "$ENV_FILE" ]]; then
         if [[ -f "workflow/envs/DiscoVir.yaml" ]]; then
            ENV_FILE="workflow/envs/DiscoVir.yaml"
         else
            echo -e "${red}[ERROR]${nc} Environment file $ENV_FILE not found!"
            exit 1
         fi
    fi

    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "       > Environment found: ${green}Yes${nc}"
        echo -e "${blu}[INFO]${nc} Activating environment from ${ylo}$ENV_FILE${nc}..."
        #run_with_spinner conda env update --name "$ENV_NAME" --file "$ENV_FILE" --prune --quiet
    else 
        echo -e "       > Environment found: ${red}No${nc}"
        echo -e "${blu}[INFO]${nc} Creating environment from ${ylo}$ENV_FILE${nc}..."
        run_with_spinner conda env create --name "$ENV_NAME" --file "$ENV_FILE" --quiet
    fi

    echo -ne "${blu}[INFO]${nc} Activating environment... "
    conda activate "$ENV_NAME"
    if [[ $? -eq 0 ]]; then 
        echo -e "${green}Active${nc}"
    else 
        echo -e "${red}Failed${nc}"
        exit 1
    fi
}

###################################
# --- Resource Management --- #
###################################

setup_resources(){
    # Diretório base de recursos do workflow
    RES_DIR="$PROJECT_DIR/resources"

    # --- 1. Diamond DB ---
    # Verifica se existe um DB externo definido
    if [[ -n "$diamond_db" ]]; then
        # Verifica se o arquivo ou prefixo existe (glob simples)
        if ls ${diamond_db}* 1> /dev/null 2>&1; then
            echo -e "${blu}[INFO]${nc} Linking external Diamond DB files to resources/diamond/..."
            
            mkdir -p "$RES_DIR/diamond"
            
            # Limpa links antigos
            rm -f "$RES_DIR/diamond/"*
            
            # Linka TODOS os arquivos que começam com o prefixo fornecido
            # Ex: se user passar /path/to/nr, vai linkar nr.00.acc, nr.00.phr, etc.
            ln -sf ${diamond_db}* "$RES_DIR/diamond/"
        else
            echo -e "${red}[ERROR]${nc} Diamond DB files not found for prefix: $diamond_db"
            exit 1
        fi
    fi

    # --- 2. Kraken2 DB ---
    # O config.yaml deve apontar para: "resources/kraken2/"
    if [[ -n "$kraken2_db" ]]; then
        # Remove a barra final se houver
        kraken2_db=${kraken2_db%/}

        if [[ -d "$kraken2_db" ]]; then
            echo -e "${blu}[INFO]${nc} Linking external Kraken2 DB to resources/kraken2/..."
            
            # Remove o diretório ou link 'kraken2' existente dentro de resources
            # Atenção: Isso substitui a pasta local pelo link para a externa
            rm -rf "$RES_DIR/kraken2_db"
            
            # Cria o link simbólico do diretório inteiro
            ln -sfn "$kraken2_db" "$RES_DIR/kraken2"
        else
            echo -e "${red}[ERROR]${nc} Kraken2 directory not found: $kraken2_db"
            exit 1
        fi


}


###################################
# --- Sample Generation --- #
###################################

generate_sample_list(){
    echo -ne "${blu}[INFO]${nc} Generating sample list..."
    sample_list="$output/samples_detected.txt"
    mkdir -p "$output"
    
    if [[ -d "$input" ]]; then
        find "$input" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) \
        | sed 's|.*/||' \
        | sed -E 's/(_1|_2|_R1|_R2|_unp|_orphans)?(_001)?\.(fastq\.gz|fastq|fasta|fa|fas)$//' \
        > "$sample_list"
    else
        touch "$sample_list"
    fi

    if [[ -n "$sra" && -f "$sra" ]]; then
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions..."
        echo "" >> "$sample_list" 
        cat "$sra" >> "$sample_list"
        echo "" >> "$sample_list"
    fi

    sort -u "$sample_list" | sed '/^$/d' > "${sample_list}.tmp" && mv "${sample_list}.tmp" "$sample_list"
    count=$(wc -l < "$sample_list")
    
    if [[ "$count" -eq 0 ]]; then
        echo -e " ${red}FAILED${nc}"
        echo -e "${red}Error: No samples found locally and no SRA file provided.${nc}"
        exit 1
    else
        echo -e " ${green}OK ($count unique samples)${nc}"
    fi
}

###########################
# --- Main Execution --- #
###########################
# Whole path of WORDDIR
PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
cd "$PROJECT_DIR"
set -e -o pipefail

# Defaults
input="data"
jobs=15
profile="profile_slurm"
temp_dir=""
snakemake_workdir=""    # Set via --workdir; defaults to PROJECT_DIR after arg parsing

# Module Defaults (Must match the keys in your config.yaml)
mod_keep_download="true"   # lowercase for yaml
mod_assembly="false"
mod_kraken2="false"
mod_diamond="false"
mod_darkmatter="false"
mod_rvdb="false"
mod_reads_kraken2="false"
mod_reads_diamond="false"
mod_download_only="false"

# --- Argument Parsing --- #
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --sra) sra="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        --profile) profile="$2"; shift 2 ;;
        
        ## Database Overrides
        --diamond_db) diamond_db="$2"; shift 2 ;;
        --kraken2_db) kraken2_db="$2"; shift 2 ;;

        ## Temp dir
        --temp-dir) temp_dir="$2"; shift 2 ;;

        ## Snakemake working directory (where SLURM jobs run from)
        --workdir) snakemake_workdir="$2"; shift 2 ;;

        # --- Module Toggles ---
        --assembly) mod_assembly="true"; shift ;;
        --kraken2) mod_kraken2="true"; shift ;;
        --diamond) mod_diamond="true"; shift ;;
        --darkmatter) mod_darkmatter="true"; shift ;;
        --rvdb) mod_rvdb="true"; shift ;;
        --download-only) mod_download_only="true" ; mod_keep_download="true"; shift ;;

        # --- Negative Flags ---
        --remove-download) mod_keep_download="false"; shift ;;
        --reads-kraken) mod_reads_kraken2="true"; shift ;;
        --reads-diamond) mod_reads_diamond="true"; shift ;;

        # --- Others --- #
        -h|--help) help; exit 0 ;;
        -v|--version) version; exit 0 ;;
        *) echo -e "${red}[ERROR]${nc} Unknown argument: $1"; help; exit 1 ;;
    esac
done

if [[ -z "$output" ]]; then 
    echo -e "${red}[ERROR]${nc} Output directory is required. Use --output <DIR>"
    exit 1
fi

manage_environment
setup_resources
generate_sample_list

echo -e "${blu}[INFO]${nc} Initializing DiscoVir Workflow..."

# 1. Locate the STATIC Main Config
main_config="$PROJECT_DIR/config/config.yaml"
if [[ ! -f "$main_config" ]]; then
    # Fallback check if it's in root
    if [[ -f "config.yaml" ]]; then
        main_config="config.yaml"
    else
        echo -e "${red}[ERROR]${nc} Could not find config.yaml in config/ or current dir."
        exit 1
    fi
fi

# 2. Generate the DYNAMIC Override Config
# This file only contains what changes per run (inputs, outputs, modules)

# Discover the true tmpdir from the selected profile if the user didn't explicitly pass --temp-dir
if [[ -z "$temp_dir" ]]; then
    echo -e "${blu}[INFO]${nc} Discovering tmpdir from profile: ${ylo}$profile${nc}"
    temp_dir=$(python3 -c "
import yaml, sys, os
prof = '$profile'
paths = [
    os.path.join(prof, 'config.yaml'),
    os.path.join('profiles', prof, 'config.yaml')
]
for p in paths:
    if os.path.exists(p):
        try:
            with open(p) as f:
                d = yaml.safe_load(f)
                t = d.get('default-resources', {}).get('tmpdir', '')
                if t:
                    print(str(t).strip())
                    sys.exit(0)
        except Exception: pass
print('/tmp/${USER}/') # fallback
")
    # Limpa possíveis quebras de linha que possam quebrar o mkdir no bash
    temp_dir=$(echo "$temp_dir" | tr -d '\r\n')
fi

# Garante que o usuário possua sua própria pasta para evitar colisão de permissões
# caso o parâmetro recuperado (ex: /scr) não possua identificação.
if [[ "$temp_dir" != *"${USER}"* ]]; then
    temp_dir="${temp_dir%/}/${USER}"
fi

# Resolve Snakemake working directory.
# Por padrão usa o mesmo scratch que --temp-dir — ambos apontam para o mesmo lugar.
# Use --workdir separadamente só se precisar de localizações distintas.
if [[ -z "$snakemake_workdir" ]]; then
    snakemake_workdir="$temp_dir"
fi
mkdir -p "$snakemake_workdir"

run_overrides="$PROJECT_DIR/run_overrides.yaml"
timestamp=$(date +%Y%m%d_%H%M%S)
backup_overrides="$output/run_overrides_${timestamp}.yaml"

echo -e "${blu}[INFO]${nc} Generating run overrides: ${ylo}$run_overrides${nc}"
echo -e "${blu}[INFO]${nc} A copy will be saved at: ${ylo}$backup_overrides${nc}"
echo -e "${blu}[INFO]${nc} Snakemake workdir   set to: ${ylo}$snakemake_workdir${nc}"

cat <<EOF > "$run_overrides"
# Dynamic Overrides generated by DiscoVir.sh
workdir: "$snakemake_workdir"
output_dir: "$output"
data_dir: "$input"
sample_list: "$sample_list"

default-resources:
  tmpdir: "$temp_dir"
  
modules:
  keep_download: $mod_keep_download
  download_only: $mod_download_only
  assembly: $mod_assembly
  kraken2: $mod_kraken2
  diamond: $mod_diamond
  darkmatter: $mod_darkmatter
  rvdb: $mod_rvdb
  reads_kraken2: $mod_reads_kraken2
  reads_diamond: $mod_reads_diamond
EOF

cp "$run_overrides" "$backup_overrides"

# --- Workflow Execution Steps ---
# We pass BOTH config files. Snakemake loads them in order.
# The second file (overrides) updates the values of the first (main).
# guarantee it works
mkdir -p "$temp_dir"
export TMPDIR="$temp_dir"
export TEMP="$temp_dir"
export TMP="$temp_dir"

echo -e "${blu}[INFO]${nc} Global TMPDIR set to: ${ylo}$TMPDIR${nc}"

# --- Paths used by all snakemake invocations ---
JOB_ID=${SLURM_JOB_ID:-local}
SHADOW_DIR="${temp_dir}/discovir_shadow/${JOB_ID}"
mkdir -p "$SHADOW_DIR"
CONDA_DIR="$PROJECT_DIR/.snakemake/conda"
SNAKEFILE="$PROJECT_DIR/workflow/Snakefile"

echo -e "${blu}[INFO]${nc} Shadow directory  set to: ${ylo}$SHADOW_DIR${nc}"
echo -e "${blu}[INFO]${nc} Tool temp dir     set to: ${ylo}$temp_dir${nc}"
echo -e "${blu}[INFO]${nc} Scratch workdir   set to: ${ylo}$snakemake_workdir${nc}"
echo -e "${blu}[INFO]${nc} Profile           set to: ${ylo}$profile${nc}"
echo -e "${blu}[INFO]${nc} Snakefile         set to: ${ylo}$SNAKEFILE${nc}"

# configfile: no Snakefile usa caminho absoluto (via workflow.snakefile), por isso
# --directory pode ser passado sem quebrar a resolução do config.
# Os jobs SLURM rodarão a partir do scratch (snakemake_workdir) em vez do home/NFS.

echo -e "\n${green}> Snakemake: Unlocking working directory...${nc}"
snakemake --profile "$profile" --snakefile "$SNAKEFILE" \
    --directory "$snakemake_workdir" \
    --configfile "$main_config" "$run_overrides" --unlock --quiet

echo -e "\n${green}> Snakemake: Creating conda environments (if needed)...${nc}"
snakemake --profile "$profile" --snakefile "$SNAKEFILE" \
    --directory "$snakemake_workdir" \
    --configfile "$main_config" "$run_overrides" --use-conda --conda-create-envs-only --quiet

echo -e "\n${green}> Snakemake: Performing a dry-run...${nc}"
snakemake --profile "$profile" --snakefile "$SNAKEFILE" \
    --directory "$snakemake_workdir" \
    --jobs $jobs --use-conda --configfile "$main_config" "$run_overrides" --dry-run
echo "---------------------------------------------------"

echo -e "\n${green}> Snakemake: Starting main execution...${nc}"

snakemake --profile "$profile" \
    --snakefile "$SNAKEFILE" \
    --directory "$snakemake_workdir" \
    --jobs $jobs \
    --use-conda \
    --conda-prefix "$CONDA_DIR" \
    --configfile "$main_config" "$run_overrides" \
    --shadow-prefix "$SHADOW_DIR" \
    --keep-going

echo -e "\n${green}> Snakemake: Creating DAG & Report...${nc}"

snakemake --report --profile "$profile" \
    --snakefile "$SNAKEFILE" \
    --directory "$snakemake_workdir" \
    --jobs $jobs \
    --use-conda \
    --configfile "$main_config" "$run_overrides" \
    --keep-going

echo -e "\n${blu}[INFO]${nc} Deactivating ${ylo}${ENV_NAME}${nc} conda environment."
conda deactivate
echo -e "${green}✔ Done.${nc}"

echo -e "\n${green} Thank you for using DiscoVir${nc}"

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
#                              version: 12.2025                                   #
###################################################################################

# --- Color Palettes ---
blue="\033[1;34m"  
green="\033[1;32m" 
red="\033[1;31m"   
ylo="\033[1;33m"   
nc="\033[0m"       

# Nome do ambiente Conda esperado (baseado no seu arquivo workflow/envs/DiscoVir.yaml)
ENV_NAME="DiscoVir"
ENV_FILE="workflow/envs/DiscoVir.yaml"

#####################
# --- Functions --- #
#####################

# --- Spinner --- #
run_with_spinner() {
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null
    local spinner=("A" "T" "G" "C")
    local i=0
    while kill -0 $pid 2>/dev/null; do
        printf "\r\033[K [${green}${spinner[i]}${nc}] Processing...    "
        i=$(( (i+1) % 4 ))
        sleep 0.2
    done
    wait $pid
    local exit_code=$?
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
 ${blue}DiscoVir: Viral Metagenomics & 'Dusk Matter' Discovery${nc}
 Author: MSc. Matheus Cosentino | Version: 12.2025 

 Usage: 
   srun DiscoVir.sh --input <DIR> --output <DIR> [OPTIONS]

 Required Arguments:
   --input <DIR>        Directory containing raw reads (.fastq.gz) or contigs (.fasta)
   --output <DIR>       Directory where results will be saved

 Optional Arguments:
   --sra <FILE>         Text file containing SRA Accession IDs for download.
   --jobs <INT>         Number of jobs (default: 15)
   --profile <STR>      Snakemake profile (default: profile_slurm)

 Database Overrides (If you want to use external DBs instead of 'resources/'):
   --diamond_db <FILE>  External Diamond database (.dmnd)
   --kraken2 <DIR>      External Kraken2 database directory
 
 Flags:
   -h, --help           Show this help message
   -v, --version        Show version
 "
}

version(){
    echo "DiscoVir v.2025.12"
}

###################################
# --- Environment Management --- #
###################################

manage_environment(){
    echo -ne "${blue}[INFO]${nc} Checking Conda environment '${ylo}${ENV_NAME}${nc}'... "

    # Tenta localizar o conda.sh para permitir ativação dentro do script
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"

    # Verifica se o ambiente existe
    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "${green}Found${nc}"
    else
        echo -e "${red}Not Found${nc}"
        echo -e "${ylo}[INFO]${nc} Creating environment from ${ENV_FILE}. This may take a while..."
        
        if [[ -f "$workdir/$ENV_FILE" ]]; then
            run_with_spinner conda env create --file "$workdir/$ENV_FILE" --name "$ENV_NAME" --quiet
        else
            echo -e "${red}[ERROR]${nc} Environment file $workdir/$ENV_FILE not found!"
            exit 1
        fi
    fi

    # Ativa o ambiente
    echo -ne "${blue}[INFO]${nc} Activating environment... "
    conda activate "$ENV_NAME"
    if [[ $? -eq 0 ]]; then
        echo -e "${green}Active${nc}"
    else
        echo -e "${red}Failed to activate${nc}"
        exit 1
    fi
}

###################################
# --- Resource Management --- #
###################################

setup_resources(){
    # Se o usuário passou bancos externos, criamos links em resources/links/
    # e informamos ao Snakemake para usar esses links.
    
    LINK_DIR="$workdir/resources/links"
    mkdir -p "$LINK_DIR"
    
    # 1. Diamond DB Override
    if [[ -n "$diamond_db" ]]; then
        if [[ -f "$diamond_db" ]]; then
            echo -e "${blue}[INFO]${nc} Linking external Diamond DB..."
            ln -sf "$diamond_db" "$LINK_DIR/external_diamond.dmnd"
            diamond_path="$LINK_DIR/external_diamond.dmnd"
        else
            echo -e "${red}[WARNING]${nc} Provided Diamond DB not found: $diamond_db. Using default."
        fi
    fi

    # 2. Kraken2 DB Override
    if [[ -n "$kraken2" ]]; then
        if [[ -d "$kraken2" ]]; then
             echo -e "${blue}[INFO]${nc} Linking external Kraken2 DB..."
            # Remove link antigo se existir para evitar aninhamento incorreto
            rm -f "$LINK_DIR/external_kraken2"
            ln -sfn "$kraken2" "$LINK_DIR/external_kraken2"
            kraken2_path="$LINK_DIR/external_kraken2"
        else
            echo -e "${red}[WARNING]${nc} Provided Kraken2 DB directory not found: $kraken2. Using default."
        fi
    fi
}

###################################
# --- Sample Generation --- #
###################################

generate_sample_list(){
    echo -ne "${blue}[INFO]${nc} Generating sample list..."
    
    sample_list="$output/samples_detected.txt"
    mkdir -p "$output"
    
    # 1. Local Files
    if [[ -d "$input" ]]; then
        find "$input" -type f \( -name "*.fastq.gz" -o -name "*.fastq" -o -name "*.fasta" -o -name "*.fa" -o -name "*.fas" \) \
        | sed 's|.*/||' \
        | sed -E 's/(_1|_2|_R1|_R2|_unp|_orphans)?\.(fastq\.gz|fastq|fasta|fa|fas)$//' \
        | sort | uniq > "$sample_list"
    else
        touch "$sample_list"
    fi

    # 2. SRA Injection
    if [[ -n "$sra" && -f "$sra" ]]; then
        echo -e "\n${ylo}[INFO]${nc} Appending SRA Accessions..."
        cat "$sra" >> "$sample_list"
    fi

    sed -i '/^$/d' "$sample_list"
    count=$(wc -l < "$sample_list")
    
    if [[ "$count" -eq 0 ]]; then
        echo -e " ${red}FAILED${nc}"
        echo -e "${red}Error: No samples found locally and no SRA file provided.${nc}"
        exit 1
    else
        echo -e " ${green}OK ($count samples)${nc}"
    fi
}

###########################
# --- Main Execution --- #
###########################

workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)

# Parse Arguments
while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
   fi
   shift
done

if [[ -z "$output" ]]; then help; exit 1; fi
input=${input:-"data"} 
profile=${profile:-profile_slurm}
jobs=${jobs:-15}

# 1. Configura Ambiente (Cria/Ativa)
manage_environment

# 2. Prepara Recursos (Links externos se necessário)
setup_resources

# 3. Detecta Amostras
generate_sample_list

echo -e "${blue}[INFO]${nc} Initializing DiscoVir Workflow..."

# Configuração Base
config_args="data_dir=$input output_dir=$output sample_list=$sample_list"

# Overrides de Recursos (apontando para os links criados)
if [[ -n "$diamond_path" ]]; then 
    config_args="$config_args resources={diamond:\"$diamond_path\"}" 
fi
if [[ -n "$kraken2_path" ]]; then 
    config_args="$config_args resources={kraken2:\"$kraken2_path\"}" 
fi

# Comando Final
cmd="snakemake --profile $profile \
    --jobs $jobs \
    --use-conda \
    --config $config_args"

echo -e "${ylo}Running command:${nc}"
echo "$cmd"
echo "---------------------------------------------------"

eval $cmd
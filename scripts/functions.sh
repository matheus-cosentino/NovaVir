
#get work directory basis
workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd) # Get working directory

# collor palletes to be used
blue="\033[1;34m"  # blue
green="\033[1;32m" # green
red="\033[1;31m"   # red
ylo="\033[1;33m"   # yellow
nc="\033[0m"       # no color

#sppinner funcion
run_with_spinner() {
    # --- Setup and Process Start ---
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null

    # Sequence is long and flat
    local spinner=("A" "T" "G" "T" "G" "T" "T" "C" "T" "G" "A" "C" "A" "A" "C" "A" "C" "G" "A" "T" "C" "A" "A" "C" "A" "T" "G")
    local seq_len=${#spinner[@]}
    local window_size=100  

    # Associate each nucleotide with its color code
    declare -A colors=( ["A"]="$blue" ["T"]="$green" ["G"]="$ylo" ["C"]="$red" )
    
    local i=0  # Index to start the window

    # --- Spinner Loop ---
    while kill -0 $pid 2>/dev/null; do
        local colored_output=""
        
        # Determine the sequence slice (window_size nucleotides starting from index i)
        # We use modular arithmetic to wrap the sequence back to the beginning
        for (( j=0; j<window_size; j++ )); do
            local index=$(( (i + j) % seq_len ))
            local nucleotide="${spinner[index]}"
            local color="${colors[$nucleotide]}"
            
            # Append Color Code + Nucleotide + No Color Reset
            colored_output+="${color}${nucleotide}${nc}"
        done
        
        # Clear the line and print the fully colored sequence window
        printf "\r\033[K [${colored_output}] Please wait    "

        # Advance the starting index (i) by one, wrapping if necessary
        i=$(( (i+1) % seq_len ))
        sleep 0.2  # Reduced sleep for a smoother scroll effect
    done
    
    # --- Cleanup and Exit ---
    wait $pid
    local exit_code=$?
    
    printf "\r\033[K"
    if [ $exit_code -eq 0 ]; then
        echo "${green}✔ Job done!${nc}"
    else
        # Use $nc to reset color after the '✖' if it wasn't done in the original code
        echo "${red}✖ Job failed with exit code $exit_code.${nc}"
    fi
}

check_connection(){
 if ping -c 1 -W 5 google.com > /dev/null 2>&1 || \
   ping -c 1 -W 5 cloudflare.com > /dev/null 2>&1
 then
    network=0
 else
    network=1
 fi    
}

check_diamond_db(){
}

check_kraken_db(){
}


download_kraken_db(){
}

download_blast_db(){
}

build_diamond_db(){
}

check_conda(){

}
# scripts/palm_annot_run.py

import subprocess
import os
import sys

s = snakemake
palm_dir_raw = s.params.palm_annot_dir
# Garante que seja string mesmo se vier como lista do config
palm_dir = str(palm_dir_raw[0]) if isinstance(palm_dir_raw, list) else str(palm_dir_raw)
log_file = s.log[0]

# --- 1. Determine the Correct Python Executable ---
# Use the python that is currently executing this script (which should be the Conda one)
python_exe = sys.executable 

# The script to run is the external PALM script, which we run *using* the Conda Python
palm_script = os.path.join(palm_dir, "py", "palm_annot.py") 

# --- 2. Construct the Command List ---
# The command starts with the explicit Python executable
command = [
    python_exe,          # Use the explicit Conda Python executable path
    palm_script,         # The external PALM script
    "--input", s.input.fasta,
    "--seqtype", s.params.seqtype,
    "--fev", s.output.fev,
    "--rdrp", s.output.rdrp,
    "--minscore", str(s.params.minscore),
    "--threads", str(s.resources.threads),
    "--minpssmscore", str(s.params.minpssmscore)
]

<<<<<<< HEAD
# --- 3. Execute the Command ---
=======
# --- 3. Set the Custom PATH Environment Variable ---
# This is still necessary for any *other* executables PALM may call
new_path = f"{palm_dir}/bin:{palm_dir}/py:{os.environ.get('PATH', '')}"
env_vars = os.environ.copy()
env_vars['PATH'] = new_path


# --- 4. Execute the Command ---
>>>>>>> b490abfa4c6c42f3c0bcec2d6f21b0f0bc0fc077
with open(log_file, "w") as log_handle:
    try:
        subprocess.run(
            command,
            check=True,
            stderr=log_handle
        )
    except subprocess.CalledProcessError as e:
        print(f"PALM annotation failed with exit code: {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The Python executable ({python_exe}) or PALM script was not found.", file=sys.stderr)
        sys.exit(1)
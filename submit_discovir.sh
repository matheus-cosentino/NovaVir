#!/bin/bash
# submit_discovir.sh
# A user-friendly launcher that runs DiscoVir in the background on the Login Node.

# 1. Generate a unique log file name based on date/time
LOG_FILE="run_discovir_$(date +%Y%m%d_%H%M%S).log"

# 2. Print a helpful starting message
echo "=================================================="
echo "      🚀 Launching DiscoVir Orchestrator 🚀       "
echo "=================================================="
echo "Command: bash DiscoVir.sh $@"
echo "Log File: $LOG_FILE"
echo "--------------------------------------------------"

# 3. Run the pipeline using nohup (The Magic Step)
# "$@" passes all arguments (like --input, --output) directly to DiscoVir.sh
nohup bash DiscoVir.sh "$@" > "$LOG_FILE" 2>&1 &

# 4. Save the Process ID (PID) to verify it started
PID=$!
echo "✅ Pipeline is running in the background (PID: $PID)."
echo "You can close this terminal or disconnect safely."
echo ""
echo "👉 To monitor progress, run:"
echo "   tail -f $LOG_FILE"
echo ""
echo "🛑 To STOP the pipeline (Emergency Kill):"
echo "   1. Kill the orchestrator: kill $PID"
echo "   2. Cancel cluster jobs:   scancel -u $USER"
echo "=================================================="
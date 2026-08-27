

#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Get the number of repeats from input argument
if [ -z "$1" ]; then
    echo "Usage: $0 <number_of_repeats> [growthRateMultiplier1] [num_ca] [num_pa] [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]"
    exit 1
fi
NUM_REPEATS=$1
GROWTH_RATE_MULTIPLIER1=${2:-0.6}
NUM_CA=${3:-1}
NUM_PA=${4:-1}
PARTITION=${5:-long}
CPUS=${6:-16}
MEM=${7:-32GB}
WALLTIME=${8:-168:00:00}
BASE_DIR=${9:-/Candida-project-main/build/Main}
OUTPUT_ROOT=${10:-/data_production/VERTICAL_ORI}
LOG_DIR=${11:-/Candida-project-main/log}

# Loop to submit jobs for each repeat
for i in $(seq 0 $((NUM_REPEATS - 1)))
do
    bash "${SCRIPT_DIR}/one_job.sh" "$i" "$GROWTH_RATE_MULTIPLIER1" "$NUM_CA" "$NUM_PA" \
        "$PARTITION" "$CPUS" "$MEM" "$WALLTIME" "$BASE_DIR" "$OUTPUT_ROOT" "$LOG_DIR"
    echo "Submitted repeat $i"
done

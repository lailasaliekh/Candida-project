#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

usage() {
    echo "Usage: $0 <number_of_repeats> [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]"
    echo "Sweeps growthRateMultiplier1 from 0 to 1 in steps of 0.1 with 50 Candida and 50 PA cells."
}

if [ -z "${1:-}" ]; then
    usage
    exit 1
fi

NUM_REPEATS=$1
PARTITION=${2:-long}
CPUS=${3:-16}
MEM=${4:-32GB}
WALLTIME=${5:-168:00:00}
BASE_DIR=${6:-/Candida-project-main/build/Main}
OUTPUT_ROOT=${7:-/data_production/growth_rate_sweep}
LOG_DIR=${8:-/Candida-project-main/log}

NUM_CA=50
NUM_PA=50
GROWTH_RATE_VALUES=("0" "0.05" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0")

for growth_rate_multiplier in "${GROWTH_RATE_VALUES[@]}"
do
    multiplier_label="${growth_rate_multiplier/./p}"
    sweep_output_root="${OUTPUT_ROOT}/growthRateMultiplier1_${multiplier_label}"
    sweep_log_dir="${LOG_DIR}/growthRateMultiplier1_${multiplier_label}"

    mkdir -p "${sweep_log_dir}"

    echo "Submitting growthRateMultiplier1=${growth_rate_multiplier} with ${NUM_CA} Candida and ${NUM_PA} PA cells"
    bash "${SCRIPT_DIR}/submit_repeats.sh" \
        "${NUM_REPEATS}" \
        "${growth_rate_multiplier}" \
        "${NUM_CA}" \
        "${NUM_PA}" \
        "${PARTITION}" \
        "${CPUS}" \
        "${MEM}" \
        "${WALLTIME}" \
        "${BASE_DIR}" \
        "${sweep_output_root}" \
        "${sweep_log_dir}"
done

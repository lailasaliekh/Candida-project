#!/bin/bash
#SBATCH --job-name=growth
#SBATCH --output=/storage/datastore-personal/s2507701/Confined_Geometries_Project/Growth_experiment_simulations/Candida-project-rep/log/growth%j.txt
#SBATCH --error=/storage/datastore-personal/s2507701/Confined_Geometries_Project/Growth_experiment_simulations/Candida-project-rep/log/growth%j.txt
#SBATCH --partition=long
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

export OMP_NUM_THREADS="16"
ulimit -c unlimited # this is for debugging if needed

usage() {
    echo "Usage: $0 <repeat_index> [growthRateMultiplier1] [num_ca] [num_pa] [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]"
}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "$0")"

RUN_MODE=0
if [ "${1:-}" = "--run-job" ]; then
    RUN_MODE=1
    shift
fi

if [ -z "${1:-}" ]; then
    usage
    exit 1
fi

REPEAT_INDEX=$1
GROWTH_RATE_MULTIPLIER1=${2:-0.6}
NUM_CA=${3:-1}
NUM_PA=${4:-1}
PARTITION=${5:-long}
CPUS=${6:-16}
MEM=${7:-32GB}
WALLTIME=${8:-168:00:00}
BASE_DIR=${9:-/storage/datastore-personal/s2507701/Confined_Geometries_Project/Growth_experiment_simulations/Candida-project-rep/build/Main}
OUTPUT_ROOT=${10:-/storage/datastore-personal/s2507701/Confined_Geometries_Project/Growth_experiment_simulations/Candida-project-rep/data_production}
LOG_DIR=${11:-/storage/datastore-personal/s2507701/Confined_Geometries_Project/Growth_experiment_simulations/Candida-project-rep/log}
MULTIPLIER_LABEL="${GROWTH_RATE_MULTIPLIER1/./p}"
GROWTH_DIR_NAME="growthRateMultiplier1_${MULTIPLIER_LABEL}"

case "${OUTPUT_ROOT}" in
    */${GROWTH_DIR_NAME})
        EFFECTIVE_OUTPUT_ROOT="${OUTPUT_ROOT}"
        ;;
    *)
        EFFECTIVE_OUTPUT_ROOT="${OUTPUT_ROOT}/${GROWTH_DIR_NAME}"
        ;;
esac

case "${LOG_DIR}" in
    */${GROWTH_DIR_NAME})
        EFFECTIVE_LOG_DIR="${LOG_DIR}"
        ;;
    *)
        EFFECTIVE_LOG_DIR="${LOG_DIR}/${GROWTH_DIR_NAME}"
        ;;
esac

if [ "$RUN_MODE" -eq 0 ]; then
    JOB_NAME="gr${MULTIPLIER_LABEL}_Ch${NUM_CA}_unch${NUM_PA}_r${REPEAT_INDEX}"
    LOG_FILE="${EFFECTIVE_LOG_DIR}/${JOB_NAME}_%j.txt"
    mkdir -p "${EFFECTIVE_LOG_DIR}"
    sbatch \
        --job-name="${JOB_NAME}" \
        --output="${LOG_FILE}" \
        --error="${LOG_FILE}" \
        --partition="${PARTITION}" \
        --time="${WALLTIME}" \
        --cpus-per-task="${CPUS}" \
        --mem="${MEM}" \
        "${SCRIPT_PATH}" --run-job "$REPEAT_INDEX" "$GROWTH_RATE_MULTIPLIER1" "$NUM_CA" "$NUM_PA" \
        "$PARTITION" "$CPUS" "$MEM" "$WALLTIME" "$BASE_DIR" "$OUTPUT_ROOT" "$LOG_DIR"
    exit $?
fi

export OMP_NUM_THREADS="${CPUS}"

REPEAT_DIR="${EFFECTIVE_OUTPUT_ROOT}/Ch${NUM_CA}_unch${NUM_PA}/repeat${REPEAT_INDEX}/"
INPUT_FILE="${REPEAT_DIR}/input_parameters.txt"

# Create the directory for the repeat
mkdir -p "$REPEAT_DIR"

cat > "$INPUT_FILE" <<EOF
repeat_index=${REPEAT_INDEX}
growthRateMultiplier1=${GROWTH_RATE_MULTIPLIER1}
num_ca=${NUM_CA}
num_pa=${NUM_PA}
partition=${PARTITION}
cpus=${CPUS}
mem=${MEM}
time=${WALLTIME}
base_dir=${BASE_DIR}
output_root=${OUTPUT_ROOT}
effective_output_root=${EFFECTIVE_OUTPUT_ROOT}
log_dir=${LOG_DIR}
effective_log_dir=${EFFECTIVE_LOG_DIR}
job_name=gr${MULTIPLIER_LABEL}_Ch${NUM_CA}_unch${NUM_PA}_r${REPEAT_INDEX}
repeat_dir=${REPEAT_DIR}
command=${BASE_DIR}/main.out ${REPEAT_DIR} ${GROWTH_RATE_MULTIPLIER1} ${NUM_CA} ${NUM_PA}
EOF

"${BASE_DIR}/main.out" "${REPEAT_DIR}" "${GROWTH_RATE_MULTIPLIER1}" "${NUM_CA}" "${NUM_PA}"

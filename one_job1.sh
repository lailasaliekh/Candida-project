#!/bin/bash
#SBATCH --job-name=beta100_zeta100
#SBATCH --output=/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/log/newbending%j.txt
#SBATCH --error=/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/log/newbending%j.txt
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB

export OMP_NUM_THREADS="8"
ulimit -c unlimited

# Get the repeat index from the argument
REPEAT_INDEX=$1
NUM_CA=30
NUM_PA=10

# Define directories
BASE_DIR="/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/build/Main"

REPEAT_DIR="/new_bending_runs/short_run/bending100_zeta100/CA${NUM_CA}_PA${NUM_PA}/repeat${REPEAT_INDEX}/"

# Create the directory for the repeat
mkdir -p "$REPEAT_DIR"

"${BASE_DIR}/main.out" "${REPEAT_DIR}" 1 2 0 1 1 4 2e-4 "${NUM_CA}" "${NUM_PA}"

#!/bin/bash
#SBATCH --job-name=read_from_prebake
#SBATCH --output=/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/log/read_14day%j.txt
#SBATCH --error=/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/log/read_14day%j.txt
#SBATCH --partition=long
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

export OMP_NUM_THREADS="16"
ulimit -c unlimited

# Get the repeat index from the argument
REPEAT_INDEX=$1
NUM_CA=1
NUM_PA=1

# Define directories
BASE_DIR="/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/Simulation_CA_PA_copy1/build/Main"

REPEAT_DIR="/readfromfile_test/vertical_extended/CAm_PA_repeat1_f107/repeat${REPEAT_INDEX}/"

# Create the directory for the repeat
mkdir -p "$REPEAT_DIR"

"${BASE_DIR}/main.out" "${REPEAT_DIR}" 1 2 0 1 1 4 2e-4 100 100

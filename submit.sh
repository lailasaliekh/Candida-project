#!/bin/bash

#SBATCH --job-name=2type_trap
#SBATCH --output=/storage/datastore-personal/s2507701/Wall/BiofilmDES-main_confinement_Sep_attempt/Simulation_twotypes_confined_trap/log/sbatch_log%j.txt
#SBATCH --error=/storage/datastore-personal/s2507701/Wall/BiofilmDES-main_confinement_Sep_attempt/Simulation_twotypes_confined_trap/log/sbatch_log%j.txt
#SBATCH --partition=long
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
export OMP_NUM_THREADS="8"

/storage/datastore-personal/s2507701/Wall/BiofilmDES-main_confinement_Sep_attempt/Simulation_twotypes_confined_trap/build/Main/main.out test/repeat00/ 1 1 0 1 1 5 2e-4 10

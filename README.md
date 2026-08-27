# Candida Project

A C++ simulation framework for modeling spatial organisation of Candida (yeast) and Pseudomonas(bacteria) in alveoili-like confined structures.
https://www.biorxiv.org/content/10.1101/2025.03.22.644766v1

## Table of Contents

- [About](#about)
- [Building](#building)  
- [Running Simulations](#running-simulations)    
- [Slurm Scripts](#slurm-scripts)
- [Parameters](#parameters)
- [Force Models](#force-models) 
- [Simulation Structure](#simulation-structure)
- [Analysis](#analysis)
---

## About

A Discrete Element Simulation of spatial oranisation of hyphal and non-hyphal Candida albicans growing with bacteria (pseudomonas aeruginosa). 
---

## Building
To compile the code run from the terminal, this only compiles on Linux operating system. One needs CMake (v >= 3.16.0) for compilation and assumes g++ 9.3.0 or above.
```bash

./quickMake.sh
```
Remember to recompile once changing the hard-coded parameters

## Running Simulation
From the terminal 
```bash
./build/Main/main.out <output_dir> [growthRateMultiplier1]
./build/Main/main.out <output_dir> [growthRateMultiplier1] [numTypeA] [numTypeB]
```

## Slurm Scripts
In order to submit jobs to the slurm cluster we can either use the following, 
for a single parameter set submission
``` bash
./one_job.sh <repeat_number> [growthRateMultiplier1] [num_ca] [num_pa] [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]
```
For loop submitting multiple repeats
``` bash
./submit_repeats.sh <number_of_repeats> [growthRateMultiplier1] [num_ca] [num_pa] [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]
```
For a sweep over `growthRateMultiplier1 = 0, 0.5, 1.0` with `50` Candida and `50` PA cells
``` bash
./submit_growthrate_sweep.sh <number_of_repeats> [partition] [cpus] [mem] [time] [base_dir] [output_root] [log_dir]
```
To allow these to be run as a program or script we use,
``` chmod +x <script.sh> ```

## Parameters
The initial number of cells of each type can now also be changed from the command line:
```bash
./build/Main/main.out test/repeat0/ 0.6 3 2
```
To change the type of candida from hyphal to non-hyphal (yeast-locked), we simply change the linking probability in:
Line 207 in ```main.cpp```
``` C++
            if (isTypeA) {
                auto* rod = new RodShapedBacterium{
                   x, y, 0,  //-->one can change to hard coded initial positions of the cells x, y, z (z=0)
                    angle, // Random angle
                    constants::pi * 0.5, 
                    RodShapedBacterium::mAvgGrwthRate,
                    4, // Type A property
                  --->  1, // if hyphal ca, 0 if yeast-locked ca
                    Candida::mRadius
                };
                initial_conditions.push_back(rod);
                numTypeA--;
            }
```

# Candida Project

A C++ simulation framework for modeling spatial organisation of Candida (yeast) and Pseudomonas(bacteria) in alveoili-like confined structures.

## Table of Contents

- [About](#about)
- [Force Models](#force-models) 
- [Building](#building)  
- [Running Simulations](#running-simulations)    
- [Slurm Scripts](#slurm-scripts)
- [Parameters](#parameters)
- [Project Structure](#project-structure)  

---

## About

A Discrete Element Simulation of spatial oranisation of hyphal and non-hyphal Candida albicans growing with bacteria (pseudomonas aeruginosa). 
---

## Building
To compile the code run from the terminal
```bash

./quickMake.sh
```
Remember to recompile once changing the hard-coded parameters

## Running Simulation
From the terminal 
```bash
./Main/main.out <output_dir>
```

## Slurm Scripts
In order to submit jobs to the slurm cluster we can either use the following, 
for a single parameter set submission
``` bash
sbatch one_job.sh <repeat_number>
```
For loop submitting multiple repeats
``` bash
./submit_repeats.sh <number_of_repeats>
```

## Parameters
The initial number of cells of each type can be changed from the ```main.cpp``` file
``` C++
  int numTypeA = 1;      // Number of ca
  int numTypeB = 1;      // Number of pa
```
To change the type of candida from hyphal to non-hyphal (yeast-locked), we simply change the linking probability in:
Line 207 in ```main.cpp```
``` C++
            if (isTypeA) {
                auto* rod = new RodShapedBacterium{
                    3/1.17, 75/1.17, 0,  // x, y, z (z=0)
                    constants::pi * 0.5, constants::pi * 0.5, // Random angle
                    RodShapedBacterium::mAvgGrwthRate,
                    4, // Type A property
                  --->  1, // if hyphal ca, 0 if yeast-locked ca
                    Candida::mRadius
                };
                initial_conditions.push_back(rod);
                numTypeA--;
            }
```


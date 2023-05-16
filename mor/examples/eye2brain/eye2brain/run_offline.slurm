#!/bin/bash

#SBATCH -J offline       # name of the job
#SBATCH -N 1             # number of nodes
#SBATCH -n 12            # number of cores
#SBATCH -t 06:45:00      # job duration
#SBATCH -o offline%j.out # standard output
#SBATCH -e offline%j.err # standard error

# set executable directory
export EXEC_DIR=../../../../build/mor/mor/examples/eye2brain

# run offline simulation
mpirun -np 12 $EXEC_DIR/eye2brain_offline feelpp_mor_eye2brain_3dP2app --config-file=${PWD}/eye2brain-3d.cfg --display-stats 1

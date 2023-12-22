#!/bin/bash

#SBATCH -N 1             # number of nodes
#SBATCH -n 12            # number of cores
#SBATCH -t 06:45:00      # job duration
#SBATCH -J eye2brain      # name of the job
#SBATCH -o %J_eye2brain_O_dirac.out # standard output
#SBATCH -e %J_eye2brain_O_dirac.err # standard error

# set executable directory
export EXEC_DIR=../../../../build/mor/mor/examples/eye2brain

# run offline simulation
mpirun -np 12 ${EXEC_DIR}/feelpp_mor_eye2brain_3dP2app --config-file=${PWD}/eye2brain-3d.cfg  --measure-index 1


# out = [301.2190009359564,300.822849307408,300.6845317877069,300.6907332716123,300.6848286873644,0]
# h = [0.08253314319606485,0.07282756078896943,0.07480314942197873,0.048449681989969025,0.05113700080520801,0.030347056736296795]
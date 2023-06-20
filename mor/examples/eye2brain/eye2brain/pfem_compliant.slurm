#!/bin/bash

#SBATCH -J pfem          # name of the job
#SBATCH -N 1             # number of nodes
#SBATCH -n 12            # number of cores
#SBATCH -t 06:45:00      # job duration
#SBATCH -o off_compliant_pfem_%j.out # standard output
#SBATCH -e off_compliant_pfem_%j.err # standard error

# output-index=1

# set executable directory
export EXEC_DIR=../../../../build/mor/mor/examples/eye2brain

# run offline simulation
mpirun -np 12 ${EXEC_DIR}/feelpp_mor_eye2brain_3dP2app --crb.use-predefined-test-sampling 1 --config-file eye2brain/eye2brain-3d.cfg --display-stats 1
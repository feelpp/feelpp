#!/bin/bash

#SBATCH -J pfem          # name of the job
#SBATCH -N 1             # number of nodes
#SBATCH --ntasks-per-node=12
#SBATCH --threads-per-core=1   # no hyperthreading
#SBATCH -t 06:45:00      # job duration
#SBATCH -o log/%j_off_compliant_pfem_.o.log # standard output
#SBATCH -e log/%j_off_compliant_pfem_.e.log # standard error

# output-index=1

# set executable directory
export EXEC_DIR=../../../../build/default/mor/examples/eye2brain

# run offline simulation
mpiexec -bind-to core ${EXEC_DIR}/feelpp_mor_eye2brain_3dP2app --crb.use-predefined-test-sampling 1 --config-file eye2brain/eye2brain-3d.cfg --display-stats 1
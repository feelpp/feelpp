#!/bin/bash

#SBATCH -J offline       # name of the job
#SBATCH -N 1             # number of nodes
#SBATCH -n 12            # number of cores
#SBATCH -t 06:45:00      # job duration
#SBATCH -o offline%j.out # standard output
#SBATCH -e offline%j.err # standard error

# set executable directory
export EXEC_DIR=../../../../build/mor/mor/examples/eye2brain

export DIM=3
export ORDER=2

index_array=(1 2 3 4 5 6 7 8 9 10)
output_names=(pointO pointA pointB pointB1 pointC pointD pointD1 pointF pointG cornea)

# run offline simulation
for i in "${!index_array[@]}";
do
    echo Generating basis for ${output_names[$i]}
    mpirun -np 12 ${EXEC_DIR}/feelpp_mor_eye2brain_${DIM}dP${ORDER}app --config-file=${PWD}/eye2brain-${DIM}d.cfg --measure-index ${index_array[$i]} --eye2brain.run.mode 2 --eye2brain.model-name eye2brain_${DIM}dP${ORDER}_${output_names[$i]}
    printf "\n\n"
done



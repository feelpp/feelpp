#!/bin/bash

PART_LIST="2 4 8 10 12 16 20 24"
PART_LIST="1 2 4 8 10 12 16 20 24"
#PART_LIST="2 4 8 10 12 24"

#H_SIZES="0.32 0.16 0.08 0.04 0.032 0.016 0.008 0.004 0.0032 0.0016 0.0008 0.0004"
#H_SIZES="0.32"

# coarse hsizes for space convergence analysis:
H_SIZES="0.032 0.016 0.008 0.004"
# additional tests
H_SIZES="0.0064"
# hsizes for space convergence analysis:
#H_SIZES="0.0032 0.0016 0.0008 0.0004"
# hsizes for redistanciation frequency analysis:
#H_SIZES="0.0006 0.0008 0.0011 0.0016 0.0023 0.0032 0.0045 0.0064 0.0091"

INPUT_MESH=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/zalesak/2d/domain.geo
OUTPUT_DIR=/home/u2/derhovsepian/feel/meshes/toolboxes/levelset/cases/benchmark/zalesak/2d
OUTPUT_FILENAME=domain

#HOSTNAME=`hostname`

CURRENT_DIR=${PWD}

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

module purge
module load feelpp-toolboxes/develop_gcc830_openmpi402

for HSIZE in ${H_SIZES}
do
    #INPUT_MESH=/home/u2/derhovsepian/feel/toolboxes/levelset/cases/benchmark/zalesak/2d/meshes/${HSIZE}/domain.msh
    mkdir -p h_${HSIZE}
    cd h_${HSIZE}
    printf "feelpp_mesh_partitioner --dim 2 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile ${OUTPUT_FILENAME}_${HSIZE}  --gmsh.hsize ${HSIZE} \n"
    (time feelpp_mesh_partitioner --dim 2 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile "${OUTPUT_FILENAME}_${HSIZE}" --gmsh.hsize ${HSIZE}) 2>&1 | tee partitioner.log &
    cd ${OUTPUT_DIR}
done
cd ${CURRENT_DIR}
wait

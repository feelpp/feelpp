#!/bin/bash

PART_LIST="2 4 8 10 12 16 20 24"
PART_LIST="1 2 4 8 10 12 16 20 24"
H_SIZES="0.032 0.016 0.008 0.004"
INPUT_MESH=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/zalesak/3d/domain3d.geo
OUTPUT_DIR=/home/u2/derhovsepian/feel/meshes/toolboxes/levelset/cases/benchmark/zalesak/3d
OUTPUT_FILENAME=domain3d
INIT_DIR=${PWD}

mkdir -p ${OUTPUT_DIR}

for HSIZE in ${H_SIZES}
do
    cd ${OUTPUT_DIR}
    mkdir -p ${HSIZE}
    cd ${HSIZE}
    printf "feelpp_mesh_partitioner --dim 3 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile ${OUTPUT_FILENAME}_${HSIZE}  --gmsh.hsize ${HSIZE} \n"
    (time feelpp_mesh_partitioner --dim 3 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile "${OUTPUT_FILENAME}_${HSIZE}" --gmsh.hsize ${HSIZE}) 2>&1 | tee partitioner.log &
done
cd ${INIT_DIR}
wait

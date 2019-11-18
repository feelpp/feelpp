#!/bin/bash

PART_LIST="2 4 8 10 12 16 20 24"
PART_LIST="1 2 4 6 8 10 12 16 20 24"
H_SIZES="0.32"
H_SIZES="0.0032 0.0016 0.0008 0.0004"
INPUT_MESH=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/square.geo
OUTPUT_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/zalesak/2d/meshes/
OUTPUT_FILENAME=square

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

for HSIZE in ${H_SIZES}
do
    mkdir -p ${HSIZE}
    cd ${HSIZE}
    printf "feelpp_mesh_partitioner --dim 2 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile ${OUTPUT_FILENAME}_${HSIZE}  --gmsh.hsize ${HSIZE} \n"
    (time feelpp_mesh_partitioner --dim 2 --part ${PART_LIST} --ifile ${INPUT_MESH} --ofile "${OUTPUT_FILENAME}_${HSIZE}" --gmsh.hsize ${HSIZE}) 2>&1 | tee partitioner.log &
    cd ${OUTPUT_DIR}
done
wait

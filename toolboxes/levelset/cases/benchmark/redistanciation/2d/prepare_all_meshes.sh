#!/bin/bash
PARTITIONER_BIN="feelpp_mesh_partitioner"
GMSH_BIN="gmsh"
DIM=2
NCORE=1

#PART_LIST="2 4 8 10 12 16 20 24"
PART_LIST="1 2 4 6 8 10 12 16 20 24"
#H_SIZES="0.32"
#H_SIZES="0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.00005"
#H_SIZES="0.1024 0.0512 0.0256 0.0128 0.0064 0.0032 0.0016 0.0008 0.0004 0.0002 0.0001"
H_SIZES="0.0004 0.0002 0.0001"
#H_SIZES="0.0004"
INPUT_GEO=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d/square.geo
OUTPUT_ROOT=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d/meshes
OUTPUT_FILENAME=square
EXCLUSIVE=""
EXCLUSIVE="--exclusive"
EXCLUDE="--exclude=atlas5,atlas6"

for HSIZE in ${H_SIZES}
do
    OUTPUT_DIR=${OUTPUT_ROOT}/${HSIZE}
    mkdir -p ${OUTPUT_DIR}
    MSH_FILE=${OUTPUT_DIR}/${OUTPUT_FILENAME}.msh
    CMD_FILE=${OUTPUT_DIR}/cmd.sh
    echo "#!/bin/bash" > ${CMD_FILE}
    echo "cd ${OUTPUT_DIR}" >> ${CMD_FILE}
    # meshing with gmsh 
    echo "module purge" >> ${CMD_FILE}
    echo "module load openmpi/3.1.4_gcc830" >> ${CMD_FILE}
    echo "module load gmsh/4.3.0_gcc830_openmpi402" >> ${CMD_FILE}
    echo "GMSH_VERSION=\$(gmsh --version)" >> ${CMD_FILE}
    echo "echo \"gmsh version: \${GMSH_VERSION} \"" >> ${CMD_FILE}
    echo "(time ${GMSH_BIN} -${DIM} ${INPUT_GEO} -o ${MSH_FILE} -clscale ${HSIZE}) 2>&1 | tee gmsh_\${GMSH_VERSION}.log" >> ${CMD_FILE}
    # partioning mesh with the feelpp partitioner
    echo "module purge" >> ${CMD_FILE}
    echo "module load feelpp.profile" >> ${CMD_FILE}
    for PART in ${PART_LIST}
    do
        echo "(time ${PARTITIONER_BIN} --dim ${DIM} --part ${PART} --ifile ${MSH_FILE} --ofile ${OUTPUT_FILENAME}) 2>&1 | tee partitioner_part_${PART}.log" >> ${CMD_FILE}
    done
    echo "(time ${PARTITIONER_BIN} --dim ${DIM} --part ${PART_LIST} --ifile ${MSH_FILE} --ofile ${OUTPUT_FILENAME}) 2>&1 | tee partitioner_all_parts.log" >> ${CMD_FILE}    
    chmod a+x ${CMD_FILE}
    #cat ${CMD_FILE}
    #printf "sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --job-name=\"2D_gmsh+partitioner_square_${HSIZE}_zalesak_redistanciation_benchmark\" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00 ${CMD_FILE} \n" 
    sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --job-name="2D_gmsh+partitioner_square_${HSIZE}_zalesak_redistanciation_benchmark" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00 ${CMD_FILE} 

done
wait


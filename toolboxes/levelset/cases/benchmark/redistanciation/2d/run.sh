#!/bin/bash
FM_BIN=/ssd/derhovsepian/feelpp//Release_build_benchmark_zalesak/toolboxes/levelset/cases/benchmark/redistanciation/feelpp_toolbox_benchmark_fast_marching_2d
NARGS=3
CFG_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d
LOGPATH=logs/redistanciation
EXPORTER_EXPORT=""
#EXPORTER_EXPORT="--exporter.export=0"
BIND_TO_CORE=""
#BIND_TO_CORE="--bind-to-core"# /!\ deprecated ! "
BIND_TO_CORE="--bind-to core"
mkdir -p ${LOGPATH} 

usage(){
    echo "Usage: $0 h np f"
    echo "h : the mesh size"
    echo "np : number of mpi processes to use"
    echo "f : file to write redistanciation time"
}

if [ $# -lt ${NARGS} ]; then
    echo "Too few arguments ($# instead of ${NARGS})"
    usage
    exit 1
fi
if [ $# -gt ${NARGS} ]; then
    echo "Too many arguments ($# instead of ${NARGS})"
    usage
    exit 1
fi

HSIZE=$1
NCORE=$2
TIMERFILE=$3

CMD="(time mpirun -np ${NCORE} ${BIND_TO_CORE} ${FM_BIN} \
    --config-file=${CFG_DIR}/redistanciation.cfg \
    --levelset.mesh.filename=${CFG_DIR}/meshes/${HSIZE}/square_p${NCORE}.json \
    --levelset.gmsh.hsize=${HSIZE} \
    --directory=toolboxes/levelset/benchmark/redistanciation/2d/h_${HSIZE} \
    --levelset.redistanciation-timer-file=${TIMERFILE} \
    ${EXPORTER_EXPORT} \
    ) 2>&1 | tee ${LOGPATH}/h_${HSIZE}_np_${NCORE}.log"

printf "${CMD} \n"
printf "${CMD} \n" >> sc_cmd.txt
eval "${CMD}"

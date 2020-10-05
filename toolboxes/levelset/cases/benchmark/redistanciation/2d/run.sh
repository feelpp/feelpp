#!/bin/bash
FM_BIN=/ssd/derhovsepian/feelpp//Release_build_benchmark_zalesak/toolboxes/levelset/cases/benchmark/redistanciation/feelpp_toolbox_benchmark_fast_marching_2d
NARGS=4
CFG_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d
EXPORTER_EXPORT=""
#EXPORTER_EXPORT="--exporter.export=0"
BIND_TO_CORE=""
#BIND_TO_CORE="--bind-to-core"# /!\ deprecated ! "
BIND_TO_CORE="--bind-to core"

usage(){
    echo "Usage: $0 c h np f"
    echo "c : the initial condition [slotteddisk,twoconcentriccircles,fourconcentriccircles]"
    echo "h : the mesh size"
    echo "np : the number of mpi processes to use"
    echo "f : the path to a directory to write the redistanciation times"
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

INITIAL_LS=$1
HSIZE=$2
NCORE=$3
TIMERDIR=$4

JSONFILE=${INITIAL_LS}.json
LOGPATH=logs/redistanciation/${INITIAL_LS}
TIMERPATH=${TIMERDIR}/${INITIAL_LS}
TIMERFILE=${TIMERPATH}/redistanciation_2d_times_${INITIAL_LS}_${HSIZE}.csv

mkdir -p ${TIMERPATH}
mkdir -p ${LOGPATH} 

CMD="(time mpirun -np ${NCORE} ${BIND_TO_CORE} ${FM_BIN} \
    --config-file=${CFG_DIR}/redistanciation.cfg \
    --levelset.filename=${CFG_DIR}/${JSONFILE} \
    --levelset.mesh.filename=${CFG_DIR}/meshes/${HSIZE}/square_p${NCORE}.json \
    --levelset.gmsh.hsize=${HSIZE} \
    --directory=toolboxes/levelset/benchmark/redistanciation/2d/${INITIAL_LS}/h_${HSIZE} \
    --levelset.redistanciation-timer-file=${TIMERFILE} \
    ${EXPORTER_EXPORT} \
    ) 2>&1 | tee ${LOGPATH}/${INITIAL_LS}_h_${HSIZE}_np_${NCORE}.log"

printf "${CMD} \n"
printf "${CMD} \n" >> sc_cmd.txt
eval "${CMD}"

#!/bin/bash
#LS_TOOLBOX_BIN=/ssd/derhovsepian/feelpp/build_develop/toolboxes/levelset/feelpp_toolbox_levelset_2d
LS_TOOLBOX_BIN=/ssd/derhovsepian/feelpp/install_develop_rpath_custom_ls_toolbox/bin/feelpp_toolbox_levelset_2d
#LS_TOOLBOX_BIN=/ssd/derhovsepian/feelpp/Debug_install_develop_rpath_custom_ls_toolbox/bin/feelpp_toolbox_levelset_2d
NARGS=3
NCORE=1
CFG_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/benchmark/zalesak/2d
HSIZE=0.004 #0.05
STABILIZATION_METHOD=SUPG
LOGPATH=logs/time_convergence
TIME_FINAL=628.0
#TIME_FINAL=0.5
TIME_SCHEME=BDF2
LS_REINIT_EVERY="-1"
QUAD_ORDER=1
QUAD_ORDER=2
#QUAD_ORDER=3
#QUAD_ORDER=4
#QUAD_ORDER=5
#QUAD_ORDER=6
#HN=${hostname 2>&1}
EXPORTER_EXPORT=""
#EXPORTER_EXPORT="--exporter.export=0"
mkdir -p ${LOGPATH} 

usage(){
    echo "Usage: $0 ts dt np"
    echo "ts : time scheme, either \"Euler\" or \"BDF2\""
    echo "dt : the time step value"
    echo "np : number of mpi processes to use"
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


TIME_SCHEME=$1
DELTA_T=$2
NCORE=$3

if [ ${TIME_SCHEME} = Euler ];
then
    ORDER=1
else
    if [ ${TIME_SCHEME} = BDF2 ];
    then
        ORDER=2
    else
        echo "Invalid time scheme (Euler or BDF2)"
        exit 1
    fi
fi


printf "h = $HSIZE, dt = ${DELTA_T}\n"
CMD="(time mpirun -np ${NCORE} ${LS_TOOLBOX_BIN} \
    --config-file=${CFG_DIR}/slotteddisk.cfg \
    --levelset.mesh.filename=${CFG_DIR}/meshes/${HSIZE}/domain_0_p${NCORE}.json \
    --directory=toolboxes/levelset/zalesak2d/time_convergence/${TIME_SCHEME}_reinitevery_${LS_REINIT_EVERY}/deltat_${DELTA_T}_h_${HSIZE}/quadorder_${QUAD_ORDER} \
    --ts.time-step=${DELTA_T} \
    --ts.time-final=${TIME_FINAL} \
    --levelset.reinit-every=${LS_REINIT_EVERY} \
    --levelset.bdf.order=${ORDER} \
    --levelset.gmsh.hsize=${HSIZE} \
    --levelset.stabilization.method=${STABILIZATION_METHOD} \
    --levelset.quad.order=$QUAD_ORDER \
    ${EXPORTER_EXPORT} \
    ) 2>&1 | tee ${LOGPATH}/${TIME_SCHEME}_reinitevery_${LS_REINIT_EVERY}_${STABILIZATION_METHOD}_h_${HSIZE}_deltat_${DELTA_T}_quadorder_${QUAD_ORDER}_np_${NCORE}.log"

printf "${CMD} \n"
printf "${CMD} \n" >> sc_cmd.txt
eval "${CMD}"




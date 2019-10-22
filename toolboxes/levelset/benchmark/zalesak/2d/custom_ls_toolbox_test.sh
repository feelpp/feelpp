#!/bin/bash
#LS_TOOLBOX_BIN=/ssd/derhovsepian/feelpp/build_develop/toolboxes/levelset/feelpp_toolbox_levelset_2d
LS_TOOLBOX_BIN=/ssd/derhovsepian/feelpp/install_develop_rpath_custom_ls_toolbox/bin/feelpp_toolbox_levelset_2d
NARGS=3
CFG_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/benchmark/zalesak/2d
EXPERIMENT_PREFIX=""
#EXPERIMENT_PREFIX="full_time_min_and_max_modGradPhi_"
LOGPATH=logs/${EXPERIMENT_PREFIX}space_convergence
TIME_FINAL=628.0
#TIME_FINAL=100.0
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
RESTART=""
#RESTART="--ts.restart=true --ts.restart.at-last-save=true"
EXPORTER_EXPORT=""
#EXPORTER_EXPORT="--exporter.export=0"
mkdir -p ${LOGPATH} 

usage(){
    echo "Usage: $0 sm h np"
    echo "sm : stabilization method, either \"GALS\", \"SUPG\", \"SGS\" or \"CIP\""
    echo "h : the mesh size"
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

STABILIZATION_METHOD=$1
HSIZE=$2
NCORE=$3

#DELTA_T=0.0
# We use dt = C*h/Umax
# with C=0.8 and Umax=0.007
# Umax is the maximum of
# the velocity field imposed
case ${HSIZE} in
    0.32 )
        DELTA_T=36.57
        ;;
    0.16 )
        DELTA_T=18.29
        ;;
    0.1 )
        DELTA_T=11.43
        ;;
    0.08 )
        DELTA_T=9.14
        ;;
    0.04 )
        DELTA_T=4.57
        ;;
    0.032 )
        DELTA_T=3.657
        ;;
    0.016 )
        DELTA_T=1.829
        ;;
    0.008 )
        DELTA_T=0.914
        ;;
    0.004 )
        DELTA_T=0.457
        ;;
    0.01 )
        DELTA_T=1.143
        ;;
    0.0032 )
        DELTA_T=0.3657
        ;;
    0.0016 )
        DELTA_T=0.1829
        ;;
    0.0008 )
        DELTA_T=0.0914
        ;;
    0.0004 )
        DELTA_T=0.0457
        ;;
esac

if [ ${TIME_SCHEME} = Euler ];
then
    ORDER=1
else
    if [ ${TIME_SCHEME} = BDF2 ];
    then
        ORDER=2
    else
        echo "Invalid time scheme (Euler or BDF2)"
        usage
        exit 1
    fi
fi

printf "h = $HSIZE, dt = ${DELTA_T}\n"
CMD="(time mpirun -np ${NCORE} ${LS_TOOLBOX_BIN} \
    --config-file=${CFG_DIR}/slotteddisk.cfg \
    --levelset.mesh.filename=${CFG_DIR}/meshes/${HSIZE}/domain_0_p${NCORE}.json \
    --directory=toolboxes/levelset/zalesak2d/${EXPERIMENT_PREFIX}space_convergence/${STABILIZATION_METHOD}_reinitevery_${LS_REINIT_EVERY}/h_${HSIZE}/quadorder_${QUAD_ORDER} \
    --ts.time-step=${DELTA_T} \
    --ts.time-final=${TIME_FINAL} \
    --levelset.reinit-every=${LS_REINIT_EVERY} \
    --levelset.ts.order=${ORDER} \
    --levelset.gmsh.hsize=${HSIZE} \
    --levelset.stabilization.method=${STABILIZATION_METHOD} \
    --levelset.quad.order=$QUAD_ORDER \
    ${EXPORTER_EXPORT} \
    ${RESTART} \
    ) 2>&1 | tee ${LOGPATH}/${TIME_SCHEME}_reinitevery_${LS_REINIT_EVERY}_${STABILIZATION_METHOD}_h_${HSIZE}_deltat_${DELTA_T}_quadorder_${QUAD_ORDER}_np_${NCORE}.log"

printf "${CMD} \n"
printf "${CMD} \n" >> cmd.txt
eval "${CMD}"




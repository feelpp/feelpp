LS_TOOLBOX_BIN=feelpp_toolbox_levelset_2d
NARGS=3
NCORE=1
TIME_FINAL=628.0 #8.0 #628.0
HSIZE=0.004 #0.05
STABILIZATION_METHOD=SUPG

#LS_REINIT_EVERY="-1"
LS_REINIT_EVERY=10
#LS_REINIT_EVERY=1

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
        echo "Invalid time scheme, use either Euler or BDF2)"
        exit 1
    fi
fi

mpirun -np ${NCORE} ${LS_TOOLBOX_BIN} \
    --config-file=slotteddisk.cfg \
    --directory=toolboxes/levelset/zalesak2d/time_convergence/${TIME_SCHEME}_reinitevery_${LS_REINIT_EVERY}/deltat_${DELTA_T} \
    --ts.time-step=${DELTA_T} \
    --ts.time-final=${TIME_FINAL} \
    --levelset.reinit-every=${LS_REINIT_EVERY} \
    --levelset.bdf.order=${ORDER} \
    --levelset.gmsh.hsize=${HSIZE} \
    --levelset.stabilization.method=${STABILIZATION_METHOD}



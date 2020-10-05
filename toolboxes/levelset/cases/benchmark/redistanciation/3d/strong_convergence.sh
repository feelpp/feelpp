RUN_SCRIPT=./run.sh
ROOT_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/3d

NCORE=1

# reminder of available sizes
#H_SIZES="0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.00005"
#H_SIZES="0.1024 0.0512 0.0256 0.0128 0.0064 0.0032 0.0016 0.0008 0.0004 0.0002 0.0001"

# We use only one size (strong convergence analysis)
HSIZE="0.0512"
HSIZE="0.0256"
HSIZE="0.0128"
HSIZE="0.01"
#HSIZE="0.0064"
#HSIZE="0.0032"

N_CORES="1 2 4 6 8 10 12 16 20 24"
N_CORES="24 20 16 12 10 8 6 4 2 1"
#N_CORES="24 20 16 12"
#N_CORES="1"
#N_CORES="24"
EXCLUDE=""
EXCLUDE="--exclude=atlas5"
#EXCLUDE="--exclude=atlas4,atlas5"
EXCLUSIVE=""
EXCLUSIVE="--exclusive"
for NCORE in ${N_CORES}
do
    TIMER_FILE=${ROOT_DIR}/redistanciation_3d_timer_h_${HSIZE}.csv
    CMD="${RUN_SCRIPT} ${HSIZE} ${NCORE} ${TIMER_FILE}"
    #eval "${CMD} &"
    #printf "sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --ntasks-per-core 1 --job-name=\"3D_h_${HSIZE}_np_${NCORE}_Z_redistanciation_benchmark\" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00${CMD}" 
    sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --ntasks-per-core 1 --job-name="3D_h_${HSIZE}_np_${NCORE}_Z_redistanciation_benchmark" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00 ${CMD} 

done

wait




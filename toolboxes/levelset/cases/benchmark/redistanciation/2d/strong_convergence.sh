#INITIAL_LS=slotteddisk
INITIAL_LS=oneconcentriccircles
#INITIAL_LS=twoconcentriccircles
#INITIAL_LS=fourconcentriccircles
RUN_SCRIPT=./run.sh
ROOT_DIR=/home/u2/derhovsepian/git/feelpp/toolboxes/levelset/cases/benchmark/redistanciation/2d
#NCORE=1

# reminder of available sizes
#H_SIZES="0.1 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.00005"
#H_SIZES="0.1024 0.0512 0.0256 0.0128 0.0064 0.0032 0.0016 0.0008 0.0004 0.0002 0.0001"

# We use only size (strong convergence analysis)
#HSIZE="0.0128"
#HSIZE="0.0064"
#HSIZE="0.0032"
#HSIZE="0.0016"
#HSIZE="0.0008"
HSIZE="0.0004"
#HSIZE="0.0002"

#N_CORES="1 2 4 6 8 10 12 16 20 24"
N_CORES="24 20 16 12 10 8 6 4 2 1"
#N_CORES="24 20 16 12 10 8 6"
#N_CORES="4 2 1"
#N_CORES="24"
#N_CORES="12"
EXCLUDE=""
EXCLUDE="--exclude=atlas5,atlas6"
#EXCLUDE="--exclude=atlas4,atlas5"
#EXCLUDE="--exclude=atlas1,atlas4,atlas5"
#EXCLUDE="--exclude=atlas1,atlas3,atlas4,atlas5"
EXCLUSIVE=""
EXCLUSIVE="--exclusive"
for NCORE in ${N_CORES}
do
    TIMER_DIR=${ROOT_DIR}/results
    CMD="${RUN_SCRIPT} ${INITIAL_LS} ${HSIZE} ${NCORE} ${TIMER_DIR}"
    #eval "${CMD} &"
    #printf "sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --ntasks-per-core 1 --job-name=\"2D_${INITIAL_LS}_h_${HSIZE}_np_${NCORE}_Z_redistanciation_benchmark\" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00${CMD}\n" 
    sbatch -n ${NCORE} -N 1-1 ${EXCLUSIVE} ${EXCLUDE} --ntasks-per-core 1 --job-name="2D_${INITIAL_LS}_h_${HSIZE}_np_${NCORE}_Z_redistanciation_benchmark" --mail-type=END --mail-user=derhovsepian@math.unistra.fr -t 48:00:00 ${CMD} 

done

wait




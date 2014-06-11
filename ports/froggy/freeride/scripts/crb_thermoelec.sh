#!/bin/bash
#OAR -n sectorHelix
#OAR -l /nodes=18,walltime=8:00:00
#OAR -p memcore=8
#OAR --project feelpp-freeride
#OAR --stdout sectorHelix.out
#OAR --stderr sectorHelix.err

source /applis/ciment/v2/env.bash
module load openmpi/1.4.4_gcc-4.6.2
echo "env loaded"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib64

cd $OAR_WORKDIR
# Number of cores
nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`
#Name of the first node
firstnode=`head -1 $OAR_NODE_FILE`
#Number of cores allocated on the first node (it is the same on all the nodes)
pernode=`grep "$firstnode\$" $OAR_NODE_FILE|wc -l`
echo "nbcores=" $nbcores
echo "nbnodes=" $nbnodes

cat $OAR_NODE_FILE

# launch the application
mpirun -np $nbcores -x LD_LIBRARY_PATH --prefix $openmpi_DIR --machinefile $OAR_NODE_FILE  \
 -mca plm_rsh_agent "oarsh" \
 $HOME/feelpp_build/boost-1.52/research/hifimagnet/applications/CRB/ThermoElectric/crb_thermoelectricCRB3DP1_nonlinearapp \
 --config-file=HRsector.cfg

# retreive outputs??

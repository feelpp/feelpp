#!/bin/bash
#OAR -n bitter
#OAR -l /nodes=4,walltime=0:30:00
#OAR --project feelpp-freeride
#OAR --stdout bitter-1.out
#OAR --stderr bitter-1.err

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
 --config-file=bitter-1.cfg

# retreive outputs??

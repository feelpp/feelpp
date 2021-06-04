#!/bin/bash

#SBATCH -p public
# number of cores
#SBATCH -n 24 
# Hyperthreading is enabled on irma-atlas, if you do not want to use it
# You must specify the following option
#SBATCH --ntasks-per-core 1
# min-max number of nodes
##SBATCH -N 1-2
# max time of exec (will be killed afterwards)
##SBATCH -t 02:00:00
# number of tasks per node
#SBATCH --tasks-per-node 24
# specify execution constraitns
##SBATCH --constraint \"intel\"
# min mem size
##SBATCH --mem=16684
# display info about cpu binding
##SBATCH --cpu_bind=verbose
# send a mail at the end of the exec
#SBATCH --mail-type=END
#SBATCH --mail-user=sala@math.unistra.fr

# If you want to have access to Feel++ logs
# export the FEELPP_SCRATCHDIR variable to an NFS mounted directory
export FEELPP_WORKDIR=/ssd/sala/job.$SLURM_JOB_ID
export FEELPP_SCRATCHDIR=/ssd/sala/job.$SLURM_JOB_ID/log

#################### OPTIONAL: 
# In case you want to use modules.
# You first have to activate the module command
source /etc/profile.d/modules.sh

# Source the configuration for Feel++ or your custom configuration
PREVPATH=`pwd`
cd /data/software/config/etc
source feelpprc.sh
cd ${PREVPATH}

# Load modules here
# This is an example of module to load
#module load gcc610.profile
#module load tools/CMake/3.5.2
module load feelpp.profile
#module load feelpp-toolboxes/develop_gcc640_openmpi1107
#module load feelpp-lib/develop_gcc640_openmpi1107

#################### OPTIONAL:

# Finally launch the job
# mpirun of openmpi is natively interfaced with Slurm
# No need to precise the number of processors to use
cfgdir=/ssd/sala/develop_before_smai/build/toolboxes/hdg/test3d0d/
###pcdir=/home/u2/sala/research/hemotumpp_src/config/preconditioners/

cd /ssd/sala/develop_before_smai/build/toolboxes/hdg
mpirun --bind-to core ./feelpp_toolbox_hdg_coupledpoisson --config-file $cfgdir/test3d0d.cfg --ts.time-step=0.01 --directory=test3d0d_0.01-bis/
# -mca btl tcp,self 


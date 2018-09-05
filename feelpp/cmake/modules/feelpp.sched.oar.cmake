###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#             Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2013-04-26
#
#  Copyright (C) 2013 Université de Strasbourg
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#

if (FEELPP_ENABLE_SCHED_OAR )
   string (REPLACE ";" " " _TMP_MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS}")
    if (FEELPP_MACHINE_NAME MATCHES "rheticus")
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.oar "#! /bin/bash
## Set Job name to app name
#OAR -n ${execname}
## Number of tasks to use and Elapsed time limit of the job
#OAR -l core=64, walltime=20:00:00
## Standard output
#OAR -O ${execname}_%jobid%.out
## Error output
#OAR -E ${execname}_%jobid%.err
## Launch job on cluster queue by default
#OAR -p cluster
## Uncomment to submit a job on the SMP machine
##OAR -p smp='YES' AND nodetype='SMP2Tb'
## Choose queue between development, short, medium & long
#OAR -q medium

# Source modules for cemracs
source /softs/cemracs_2015/cemracs.sh

# To display some info:
# Number of cores
#nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
#nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`
#Name of the first node
#firstnode=`head -1 $OAR_NODE_FILE`

#Number of cores allocated on the first node (it is the same on all the nodes)
#pernode=`grep \"$firstnode\$\" $OAR_NODE_FILE|wc -l`
#echo \"nbcores=\" $nbcores
#echo \"nbnodes=\" $nbnodes

# Important note:
# If you export additional variables like FEELPP_WORKDIR and FEELPP_SCRATCHDIR
# you need to export them for mpirun with -x
# We ended up having problem with Ginac when not exporting those variables
# (They were only set for for several processes)
# sample exports
#export FEELPP_WORKDIR=/home/user/feel/job.$OAR_JOBID
#export FEELPP_SCRATCHDIR=/home/user/logfiles/job.$OAR_JOBID

# launch the application
# For OpenMPI 1.6, use orte_rsh_agent instead of plm_rsh_agent
# See https://www.grid5000.fr/mediawiki/index.php/Run_MPI_On_Grid%275000
# See http://oar.imag.fr/docs/2.5/user/usecases.html
# Note: In OpenMPI 1.6, pls_rsh_agent was replaced by orte_rsh_agent. Note: In OpenMPI 1.8, orte_rsh_agent was replaced by plm_rsh_agent.

# Advice: Use absolute paths to ensure that the executables and config files are found. 
# Options that modify the paths, like --nochdir, might also be the source of errors for failing submissions
        ")

        if ( FEELPP_APP_CFG )
            foreach( cfg ${FEELPP_APP_CFG} )
                get_filename_component( CFG_NAME ${cfg} NAME )
                file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.oar 
                    "mpirun -x LD_LIBRARY_PATH -machinefile $OAR_NODEFILE ${_TMP_MPIEXEC_PREFLAGS} \\
    ${CMAKE_CURRENT_BINARY_DIR}/${execname} \\
    --config-file=${CMAKE_CURRENT_BINARY_DIR}/${cfg} "
                )
            endforeach()
        else()
            file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${execname}.oar 
            "mpirun -x LD_LIBRARY_PATH -machinefile $OAR_NODEFILE ${_TMP_MPIEXEC_PREFLAGS} \\
    ${CMAKE_CURRENT_BINARY_DIR}/${execname} "
            )
        endif()

    else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.oarsub "#! /bin/bash
#OAR -r ${execname}                # Request name
#OAR -l /nodes=64,walltime=0:30:00 # Number of tasks to use and Elapsed time limit of the job
#OAR --stdout ${execname}_%jobid%.o     # Standard output. %I is the job id
#OAR --stdout ${execname}_%jobid%.e     # Error output. %I is the job id
#OAR --project hpcfeelpp     # Project ID
##OAR --notify exec:/usr/local/bin/sendmail.sh     # Send mail upon exit

# if you use nix instead of modules comment out the 4 following lines
#set -x
source /applis/ciment/v2/env.bash
module load openmpi/1.4.4_gcc-4.6.2
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib64

# if you use nix instead of modules uncomment the next line
# source /applis/site/nix.sh
 
cd $OAR_WORKDIR
# Number of cores
nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`
#Name of the first node
firstnode=`head -1 $OAR_NODE_FILE`
#Number of cores allocated on the first node (it is the same on all the nodes)
pernode=`grep \"$firstnode\$\" $OAR_NODE_FILE|wc -l`
echo \"nbcores=\" $nbcores
echo \"nbnodes=\" $nbnodes

cat $OAR_NODE_FILE

# launch the application
mpirun -np $nbcores -x LD_LIBRARY_PATH --prefix $openmpi_DIR -machinefile $OAR_NODEFILE ${_TMP_MPIEXEC_PREFLAGS} \\
    ${CMAKE_CURRENT_BINARY_DIR}/${execname}
")
    endif()
endif()

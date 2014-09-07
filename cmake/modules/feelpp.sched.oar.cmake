###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
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
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${execname}.oarsub "#! /bin/bash
#OAR -r ${execname}                # Request name
#OAR -l /nodes=64,walltime=0:30:00 # Number of tasks to use and Elapsed time limit of the job
#OAR --stdout ${execname}_%I.o     # Standard output. %I is the job id
#OAR --stdout ${execname}_%I.e     # Error output. %I is the job id
#OAR --project feelpp-freeride     # Project ID

#set -x
source /applis/ciment/v2/env.bash
module load openmpi/1.4.4_gcc-4.6.2
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
 ${execname}
")
  endif()

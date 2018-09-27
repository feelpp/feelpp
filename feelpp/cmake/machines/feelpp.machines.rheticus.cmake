###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2015-08-05
#
#  Copyright (C) 2015 Université de Strasbourg
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

set(FEELPP_ENABLE_SCHED_OAR TRUE)
# see https://oar.readthedocs.org/en/2.5/user/usecases.html#using-mpi-with-oarsh
# option shall depend on MPI_FLAVOR (here defined only OpenMPI
set(MPIEXEC_PREFLAGS "-machinefile $OAR_NODEFILE -mca plm_rsh_agent \"oarsh\" " CACHE STRING "These flags will be directly before the executable that is being run by MPIEXEC.")
message(STATUS "on rheticus : MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} ")

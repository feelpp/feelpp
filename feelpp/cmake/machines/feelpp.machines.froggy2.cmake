###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-04-12
#
#  Copyright (C) 2013-2015 Feel++ Consortium
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
set(FEELPP_ENABLE_MANUAL OFF)
set(Boost_NO_BOOST_CMAKE TRUE)
set(Boost_NO_SYSTEM_PATHS TRUE)

set(FEELPP_ENABLE_SCHED_OAR TRUE)

# see https://oar.readthedocs.org/en/2.5/user/usecases.html#using-mpi-with-oarsh
# option shall depend on MPI_FLAVOR (here defined only OpenMPI
set(MPIEXEC_PREFLAGS -mca plm_rsh_agent \"oarsh\" CACHE STRING "These flags will be directly before the executable that is being run by MPIEXEC.")
message(STATUS "on froggy2 : MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} ")

# Work around to get glog to compile
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PERMIT_BACKWARD_HASH " )
endif()
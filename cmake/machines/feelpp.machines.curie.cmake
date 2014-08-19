###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-04-12
#
#  Copyright (C) 2013 Feel++ Consortium
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

# enable the CCC script generation
OPTION(FEELPP_ENABLE_SCHED_CCC "Enable Feel++ tgcc/ccc submission scripts generation" ON)

# disable NLOPT for now
OPTION( FEELPP_ENABLE_NLOPT "Enable NLOPT (NonLinear Optimisation Library)" OFF)

# set INTEL_ROOT to CXX_INTEL_ROOT to find MKL
set(INTEL_ROOT $ENV{CXX_INTEL_ROOT})

# find the gfortran library
FIND_LIBRARY(GFORTRAN_LIBRARY
    NAMES
    gfortran
    PATHS
    /usr/local/gcc-4.8.1/lib
    $ENV{LIBRARY_PATH}
)
message(STATUS "curie gfortran lib: ${GFORTRAN_LIBRARY} ")

# this is for petsc on curie
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS}  -L/usr/local/intel-14.0.3.174/lib/intel64/ -L/usr/local/intel-14.0.3.174/composer_xe_2013_sp1.3.174/compiler/lib/intel64 -L/usr/local/hwloc-1.7.1/lib/")
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  -L/usr/local/intel-14.0.3.174/lib/intel64/ -L/usr/local/intel-14.0.3.174/composer_xe_2013_sp1.3.174/compiler/lib/intel64 -L/usr/local/hwloc-1.7.1/lib/")
set(FEELPP_LIBRARIES "-lhwloc -lifport -lifcore -limf -lsvml -liomp5 -lirc -liomp5")

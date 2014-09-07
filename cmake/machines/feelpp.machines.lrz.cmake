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

set(FEELPP_ENABLE_SCHED_LOADLEVELER ON CACHE BOOL "")

set(FEELPP_ENABLE_MANUAL OFF)
#set(FEELPP_ENABLE_BENCHMARKS OFF)
set(FEELPP_ENABLE_OPENGL OFF)

# Disable use of standard c headers in ginac-excompiler
set(USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER OFF)

set(BLAS_blas_LIBRARY $ENV{packagesBaseDir}/blas/BLAS-gcc$ENV{gccVersion}/blas_LINUX.a)

# find the gfortran library
FIND_LIBRARY(GFORTRAN_LIBRARY
    NAMES
    gfortran
    PATHS
    /lrz/sys/compilers/gcc/$ENV{gccVersion}/lib64/
    $ENV{LIBRARY_PATH}
)
message(STATUS "lrz gfortran lib: ${GFORTRAN_LIBRARY} ")

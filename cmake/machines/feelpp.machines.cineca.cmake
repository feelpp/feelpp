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

set(FEELPP_USE_CLANG_LIBCXX ON)
set(FEELPP_USE_STATIC_LINKAGE ON)

#IF( ("${CMAKE_CXX_COMPILER_ID}" MATCHES "XL") )
set(BLAS_blas_LIBRARY /cineca/prod/libraries/blas/2007/bgq-xl--1.0/lib/libblas.a)
message(STATUS "Blas: ${BLAS_blas_LIBRARY}")
set(LAPACK_lapack_LIBRARY /cineca/prod/libraries/lapack/3.4.1/bgq-xl--1.0/lib/liblapack.a )
message(STATUS "LAPACK: ${LAPACK_lapack_LIBRARY}")
#endif()

set(MPI_CXX_COMPILER "$ENV{WORK}/local/bin/mpic++")
message(STATUS "mpi c++ compiler @cineca: ${MPI_CXX_COMPILER}")

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

set(FEELPP_USE_CLANG_LIBCXX OFF)
set(FEELPP_USE_STATIC_LINKAGE ON)
set(FEELPP_ENABLE_BUILD_STATIC ON)
#IF( ("${CMAKE_CXX_COMPILER_ID}" MATCHES "XL") )
set(BLAS_blas_LIBRARY $ENV{PETSC_DIR}/lib/libfblas.a)
message(STATUS "Blas: ${BLAS_blas_LIBRARY}")
set(LAPACK_lapack_LIBRARY $ENV{PETSC_DIR}/lib/libflapack.a)
message(STATUS "LAPACK: ${LAPACK_lapack_LIBRARY}")
#endif()
#set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS}  -fPIC -Bdynamic -dynamic -Wl,--allow-multiple-definition" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXEC_LINKER_FLAGS} -dynamic -Bdynamic -Wl,--allow-multiple-definition  -stdlib=libstdc++" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++" )

#set(MPI_CXX_COMPILER "$ENV{WORK}/local/bin/mpic++")
#message(STATUS "mpi c++ compiler @cineca: ${MPI_CXX_COMPILER}")
#set(MPI_CXX_INCLUDE_PATH "/bgsys/drivers/V1R2M1/ppc64/comm/lib/gnu;/bgsys/drivers/V1R2M1/ppc64;/bgsys/drivers/V1R2M1/ppc64/comm/sys/include;/bgsys/drivers/V1R2M1/ppc64/spi/include;/bgsys/drivers/V1R2M1/ppc64/spi/include/kernel/cnk")
#set(MPI_CXX_LIBRARIES "/bgsys/drivers/V1R2M1/ppc64/comm/lib/libmpich-gcc.so;/bgsys/drivers/V1R2M1/ppc64/comm/lib/libopa-gcc.so;/bgsys/drivers/V1R2M1/ppc64/comm/lib/libmpl-gcc.so;/bgsys/drivers/V1R2M1/ppc64/comm/lib/libpami-gcc.so;/bgsys/drivers/V1R2M1/ppc64/spi/lib/libSPI.a;/bgsys/drivers/V1R2M1/ppc64/spi/lib/libSPI_cnk.a;/gpfs/work/LI03s_DDFD/bgq/gnu-linux-4.7.2/powerpc64-bgq-linux/lib/librt.so;/gpfs/work/LI03s_DDFD/bgq/gnu-linux-4.7.2/powerpc64-bgq-linux/lib/libpthread.so;/gpfs/work/LI03s_DDFD/bgq/gnu-linux-4.7.2/powerpc64-bgq-linux/lib/libstdc++.so;/gpfs/work/LI03s_DDFD/bgq/gnu-linux-4.7.2/powerpc64-bgq-linux/lib/libpthread.so")
#set(MPI_C_LIBRARIES ${MPI_CXX_LIBRARIES})

if (FEELPP_USE_STATIC_LINKAGE )                                                                                                                                                                                                               
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)                                                                                
  set(Boost_USE_STATIC_LIBS   ON) 
endif()  


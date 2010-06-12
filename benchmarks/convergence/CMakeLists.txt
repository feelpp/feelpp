###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-12-01
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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
INCLUDE_DIRECTORIES (
  ${LIFE_SOURCE_DIR}/benchmarks/convergence
  )
ENABLE_TESTING()
foreach( BENCH laplacian1dP1 laplacian1dP5
    laplacian1d2dP1 laplacian1d2dP5
    laplacian2dP1 laplacian2dP5
    laplacian2d3dP1
    laplacian2dP1h laplacian2dP5h
    laplacian3dP1 laplacian3dP2 laplacian3dP3 laplacian3dP5
    laplacian3dP1h
    laplacian2dP1v laplacian2dP4v )

    life_add_bench( ${BENCH} )

    # test each backend
    foreach( BACKEND petsc trilinos )
      add_test( life_bench_${BENCH}-${BACKEND}-strong
        ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
        --program=${targetname}
        --path=${LIFE_BUILD_DIR}/benchmarks/convergence/
        --cmdargs="--backend=${BACKEND} --export=0 ${${BACKEND}_SOLVER_OPTIONS}" )
      add_test( life_bench_${BENCH}-${BACKEND}-weak
        ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
        --program=${targetname}
        --mode=0
        --path=${LIFE_BUILD_DIR}/benchmarks/convergence/
        --cmdargs="--backend=${BACKEND} --weak --export=0 ${${BACKEND}_SOLVER_OPTIONS}")

    endforeach()

endforeach()


foreach( BENCH
    nonlinearpow
    ale1 ale2 ale3 ale4 ale5 )

    life_add_bench( ${BENCH} )

    # test each backend
    foreach( BACKEND petsc trilinos )
      add_test( life_bench_${BENCH}-${BACKEND}
        ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
        --program=${targetname}
        --mode=0
        --path=${LIFE_BUILD_DIR}/benchmarks/convergence/
        --cmdargs="--backend=${BACKEND} --export=0 ${${BACKEND}_LU_SOLVER_OPTIONS}" )
    endforeach()

endforeach()


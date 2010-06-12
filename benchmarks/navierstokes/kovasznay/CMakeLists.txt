###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-06-12
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
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
  ${LIFE_SOURCE_DIR}/benchmarks/navierstokes/kovasznay
  )
ENABLE_TESTING()

# use ilu for the test cases
foreach( BENCH stokes11 stokes22 stokes33 stokes21 stokes31 stokes32 )

    life_add_bench( ${BENCH} )

    # test each backend
    foreach( BACKEND petsc trilinos )
      add_test( life_bench_${BENCH}-${BACKEND}
        ${LIFE_SOURCE_DIR}/benchmarks/python/batch_validation.py
        --program=${targetname}
        --mode=0
        --path=${LIFE_BUILD_DIR}/benchmarks/convergence/
        --cmdargs="--backend=${BACKEND} --export=0 ${${BACKEND}_ILU_SOLVER_OPTIONS} --penalisation=1e-3" )
    endforeach()

endforeach()


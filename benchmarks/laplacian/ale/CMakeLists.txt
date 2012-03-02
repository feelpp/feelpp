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


foreach( BENCH ale1 ale2 ale3 ale4 ale5 )

    feel_add_bench( ${BENCH} )

    # test each backend
    foreach( BACKEND ${FEELPP_BACKENDS} )
      add_test( feel_bench_${BENCH}-${BACKEND}
        ${PYTHON_EXECUTABLE} ${FEELPP_SOURCE_DIR}/benchmarks/python/batch_validation.py
        --program=${targetname}
        --mode=0
        --path=${FEELPP_BUILD_DIR}/benchmarks/laplacian/ale/
        --cmdargs="--exporter-format=gmsh --backend=${BACKEND} --export=0 ${${BACKEND}_LU_SOLVER_OPTIONS}" )
      set_property(TEST feel_bench_${BENCH}-${BACKEND}  PROPERTY LABELS benchmarks)
    endforeach()

endforeach()


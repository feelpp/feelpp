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
add_custom_target(genfiles ALL)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_custom_target(feelpp_assembly)

foreach(D 2 3)
  foreach(N 1 2 3 4)
    feelpp_add_application( assembly_${D}d_p${N} SRCS  assembly.cpp  DEFS FEELPP_DIM=${D} FEELPP_ORDER=${N})
    add_dependencies(feelpp_assembly feelpp_assembly_${D}d_p${N})
  endforeach()
endforeach()
    
add_executable( feel_perf
  bench1_main.cpp
  bench1_run.cpp
  bench1_run1d.cpp
  bench1_run2d.cpp
  bench1_run3d.cpp
  )
target_link_libraries(feel_perf feelpp ${FEELPP_LIBRARIES} )

add_executable( feel_perf2
  perf2.cpp
  )
target_link_libraries(feel_perf2 feelpp ${FEELPP_LIBRARIES} )

add_subdirectory(laplacian)
add_subdirectory(dar)
add_subdirectory(stokes)
add_subdirectory(mortar)
add_subdirectory(tilted)
add_subdirectory(curvature)

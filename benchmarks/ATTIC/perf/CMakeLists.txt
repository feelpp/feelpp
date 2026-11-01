# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#
# SPDX-FileCopyrightText: 2009 Université Joseph Fourier (Grenoble I)
# SPDX-FileCopyrightText: 2009-2026 University of Strasbourg
#
# SPDX-License-Identifier: LGPL-3.0-or-later

add_custom_target(genfiles ALL)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_custom_target(feelpp_assembly)

foreach(D 2 3)
  foreach(N 1 2 3 4)
    feelpp_add_application( assembly_${D}d_p${N} SRCS  assembly.cpp  DEFS FEELPP_DIM=${D} FEELPP_ORDER=${N} EXEC ASSEMBLY_APP)
    add_dependencies(feelpp_assembly ${ASSEMBLY_APP})
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


add_subdirectory(laplacian)
add_subdirectory(dar)
add_subdirectory(stokes)
add_subdirectory(mortar)
add_subdirectory(tilted)
add_subdirectory(curvature)
add_subdirectory(io)

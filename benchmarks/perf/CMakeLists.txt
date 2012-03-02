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
macro(OVERWITE_IF_DIFFERENT thetarget filename var dummy)
  if ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
    # be careful if file does not exist we use dummy to generate the cpp file which will
    # then be overwritten using the cmake -E copy_if_different command
    configure_file(${dummy}  ${CMAKE_CURRENT_BINARY_DIR}/${filename})
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${var})
  add_custom_command(TARGET ${thetarget} COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/copy_${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable( feel_perf
  bench1_main.cpp
  bench1_run.cpp
  bench1_run1d.cpp
  bench1_run2d.cpp
  bench1_run3d.cpp
  )
target_link_libraries(feel_perf ${FEELPP_LIBRARIES} )
# add_test( perf COMMAND feel_perf --hsi

add_executable( feel_perf2
  perf2.cpp
  )
target_link_libraries(feel_perf2 ${FEELPP_LIBRARIES} )

add_subdirectory(laplacian)
add_subdirectory(dar)
add_subdirectory(stokes)
add_subdirectory(mortar)

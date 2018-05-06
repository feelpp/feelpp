###  CMakeLists.txt; coding: utf-8 ---

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 06 Feb 2015
#
#  Copyright (C) 2015-2017 Feel++ Consortium
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
cmake_minimum_required(VERSION 3.3)

project(quickstart)
set(PROJECT_SHORTNAME "qs")

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

enable_testing()

if ( NOT TARGET check )
  add_custom_target(check)
endif()

if ( NOT TARGET quickstart )
  add_custom_target(quickstart)
endif()

#
# Laplacian
#
feelpp_add_application( laplacian_2d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_LAPLACIAN_2D INSTALL TESTS)
feelpp_add_application( laplacian_3d SRCS qs_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_LAPLACIAN_3D INSTALL TESTS)
feelpp_add_testcase( laplacian PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_LAPLACIAN_2D} ${QS_LAPLACIAN_3D}  )
add_dependencies(check ${QS_LAPLACIAN_2D} ${QS_LAPLACIAN_3D} )


#
# HDG Laplacian
#
#feelpp_add_application( hdg_laplacian_2d SRCS qs_hdg_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_LAPLACIAN_2D INSTALL TEST "--config-file" "${CMAKE_CURRENT_SOURCE_DIR}/laplacian/circle/circle-dirichlet.cfg")
#feelpp_add_application( hdg_laplacian_3d SRCS qs_hdg_laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_LAPLACIAN_3D INSTALL TEST "--config-file" "${CMAKE_CURRENT_SOURCE_DIR}/laplacian/tetrahedron/tetrahedron.cfg" "--gmsh.hsize" "4")
#add_dependencies(quickstart ${QS_HDG_LAPLACIAN_2D} ${QS_HDG_LAPLACIAN_3D}  )
#add_dependencies(check ${QS_HDG_LAPLACIAN_2D} ${QS_HDG_LAPLACIAN_3D} )

#
# ADR
#
feelpp_add_application( adr_2d SRCS qs_adr.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ADR_2D INSTALL TESTS)
feelpp_add_application( adr_3d SRCS qs_adr.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ADR_3D INSTALL NO_TEST)
feelpp_add_testcase( adr PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_ADR_2D} ${QS_ADR_3D}  )
add_dependencies(check ${QS_ADR_2D} ${QS_ADR_3D} )

#
# Stokes
#
feelpp_add_application( stokes_2d SRCS qs_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_STOKES_2D INSTALL TESTS)
feelpp_add_application( stokes_3d SRCS qs_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_STOKES_3D INSTALL TESTS)
feelpp_add_testcase( stokes PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_STOKES_2D} ${QS_STOKES_3D})
add_dependencies(check ${QS_STOKES_2D} ${QS_STOKES_3D})

#
# Linear elasticity
#
feelpp_add_application( elasticity_2d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_2D INSTALL )
feelpp_add_application( elasticity_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_3D INSTALL )
feelpp_add_testcase( elasticity PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})
add_dependencies(check ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})

add_custom_target(install-quickstart
  DEPENDS
  quickstart
  COMMAND
  "${CMAKE_COMMAND}"
  -DCMAKE_INSTALL_COMPONENT=testcases  -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Quickstart -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  -DCMAKE_INSTALL_COMPONENT=Bin -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
  )

if ( NOT TARGET install-feelpp-base )
  add_custom_target( install-feelpp-base DEPENDS install-quickstart)
  else()
  add_dependencies(install-feelpp-base install-quickstart)
endif()

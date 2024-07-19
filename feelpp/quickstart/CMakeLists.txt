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
cmake_minimum_required(VERSION 3.7.2)

project(quickstart)
set(PROJECT_SHORTNAME "qs")

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
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

if ( NOT TARGET hdg_stokes )
  add_custom_target(hdg_stokes)
endif()

add_subdirectory( python )
add_subdirectory( laplacian )

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
feelpp_add_application( hdg_stokes_2d SRCS qs_hdg_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_STOKES_2D INSTALL TESTS)
feelpp_add_application( hdg_stokes_3d SRCS qs_hdg_stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_STOKES_3D INSTALL TESTS)
feelpp_add_testcase( stokes PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_STOKES_2D} ${QS_STOKES_3D} ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D} )
add_dependencies(check ${QS_STOKES_2D} ${QS_STOKES_3D} ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D})
add_dependencies(hdg_stokes ${QS_HDG_STOKES_2D} ${QS_HDG_STOKES_3D})

#
# Linear elasticity
#
feelpp_add_application( elasticity_2d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( elasticity_3d SRCS qs_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_3D INSTALL )
feelpp_add_application( elasticity_pure_traction_2d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_ELASTICITY_PURE_TRACTION_2D INSTALL )
feelpp_add_application( elasticity_pure_traction_3d SRCS qs_elasticity_pure_traction.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_ELASTICITY_PURE_TRACTION_3D INSTALL  )
feelpp_add_application(elasticity_contact SRCS qs_elasticity_contact.cpp qs_elasticity_rigid_impl.cpp qs_elasticity_contact_static_impl.cpp qs_elasticity_contact_dynamic_impl.cpp INCLUDE_IN_ALL DEFS EXEC elasticity_contact INSTALL)


feelpp_add_testcase( elasticity PREFIX qs CATEGORY quickstart )
add_dependencies(quickstart ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})
add_dependencies(check ${QS_ELASTICITY_2D} ${QS_ELASTICITY_3D})

#
# HDG elasticity
#
feelpp_add_application( hdg_elasticity_2d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 EXEC QS_HDG_ELASTICITY_2D INSTALL TESTS )
feelpp_add_application( hdg_elasticity_3d SRCS qs_hdg_elasticity.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 EXEC QS_HDG_ELASTICITY_3D INSTALL TESTS )

add_dependencies(quickstart ${QS_HDG_ELASTICITY_2D} ${QS_HDG_ELASTICITY_3D})
add_dependencies(check ${QS_HDG_ELASTICITY_2D} ${QS_HDG_ELASTICITY_3D})

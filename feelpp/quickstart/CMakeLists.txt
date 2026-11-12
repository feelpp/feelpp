# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileContributor: Philippe Ricka <rickaphilippe@gmail.com>
#
# SPDX-FileCopyrightText: 2011-2026 University of Strasbourg
#
# SPDX-License-Identifier: GPL-2.0-or-later

cmake_minimum_required(VERSION 3.21)

project(quickstart)
set(PROJECT_SHORTNAME "qs")

if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "quickstart" )
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

enable_testing()

# Set quickstart label for all tests in this directory tree
set_directory_properties(PROPERTIES LABELS "quickstart")

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
add_subdirectory( feelsystem )
add_subdirectory( maxwell )

feelpp_get_compile_definition(Feelpp::feelpp FEELPP_HAS_PETSC_TAO)
if ( FEELPP_HAS_PETSC_TAO )
  add_subdirectory( optimization )
endif()

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

add_subdirectory( elasticity )

# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2014 Feel++ Consortium
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
cmake_minimum_required (VERSION 3.13.0)

if ( FEELPP_COMPONENT STREQUAL "toolboxes" OR 
      FEELPP_COMPONENT STREQUAL "mor" OR
      FEELPP_COMPONENT STREQUAL "testsuite" OR
      FEELPP_COMPONENT STREQUAL "benchmarks" OR
      FEELPP_COMPONENT STREQUAL "feelpp/pyfeelpp" OR
      FEELPP_COMPONENT STREQUAL "toolboxes/pyfeelpp-toolboxes" OR
      FEELPP_COMPONENT STREQUAL "mor/pyfeelpp-mor" OR
      FEELPP_COMPONENT STREQUAL "pyfeelpp-all"   )

      message(STATUS "[component] ${FEELPP_COMPONENT}")
      add_subdirectory(${FEELPP_COMPONENT})
      
else()

if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
  message(STATUS "[feelpp] use new policy CMP0022" )
  cmake_policy(SET CMP0022 NEW)
	message(STATUS "[feelpp] use NEW policy CMP0038" )
	cmake_policy(SET CMP0038 NEW)
	message(STATUS "[feelpp] use NEW policy CMP0042" )
	cmake_policy(SET CMP0042 NEW)
endif()
if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13.0 )
  message(STATUS "[feelpp] use new policy CMP0079" )
  cmake_policy(SET CMP0079 NEW)
endif()

if (POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
endif()

if (POLICY CMP0045)
# error on non-existent target in get_target_property
cmake_policy(SET CMP0045 NEW)
endif()

#
# Bug in cmake that delete destroy library path if LIBRARY_PATH is not empty
# Ugly hack that disable the environement variable
# Waiting for cmake to fix it
OPTION(FEELPP_RESET_ENV_LIBRARY_PATH "apply SET(ENV{LIBRARY_PATH} )" ON )
if (FEELPP_RESET_ENV_LIBRARY_PATH)
  SET(ENV{LIBRARY_PATH} "")
endif()

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# guard against bad build-type strings

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(     NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release"
    AND NOT cmake_build_type_tolower STREQUAL "minsizerel"
    AND NOT cmake_build_type_tolower STREQUAL "asan"
    AND NOT cmake_build_type_tolower STREQUAL "coverage"
    AND NOT cmake_build_type_tolower STREQUAL "debugfull"
    AND NOT cmake_build_type_tolower STREQUAL "none"
    AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are None, Debug, Release, MinSizeRel, RelWithDebInfo (case-insensitive), Asan, Coverage, DebugFull.")
endif()


#SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
#SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
  project (Feel++ C CXX )
else()
  project (Feel++ C CXX Fortran)
endif()


macro(set_config_option VARNAME STRING)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "[feelpp] Found " ${STRING})
endmacro(set_config_option)

# this variable controls the default value of the "FEELPP_ENABLE_XXX" options which are
# normally set to ON (useful if you want to configure a minimal version of feel++:
# e.g. "cmake -DFEELPP_DEFAULT=0 -DFEELPP_ENABLE_TESTS=1 -DFEELPP_ENABLE_=1")
set(FEELPP_DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")

OPTION(FEELPP_ENABLE_VERBOSE_CMAKE "enable Feel++ verbose cmake" OFF)
OPTION(FEELPP_ENABLE_BUILD_STATIC "Enable Feel++ library statically built" OFF )
if ( NOT FEELPP_ENABLE_BUILD_STATIC )
  set(BUILD_SHARED_LIBS ON)
endif()
OPTION(FEELPP_ENABLE_ASAN "enable asan flags" OFF)
OPTION(FEELPP_BENCHMARK_FLAGS "enable benchmarks flags" OFF)
OPTION(FEELPP_ENABLE_DOCUMENTATION "enable Feel++ documentation (aka examples)" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_BENCHMARKS "enable Feel++ benchmarks" OFF)
OPTION(FEELPP_ENABLE_TOOLBOXES "enable Feel++ toolboxes" OFF)
option(FEELPP_ENABLE_MOR "Enable MOR(Model Order Reduction) Applications in Feel++" OFF)

# enable PCH only for GCC as for now
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	OPTION(FEELPP_ENABLE_PCH_FOR_APPLICATIONS "enable Pre-Compiled Headers (PCH) for Feel++ applications" OFF)
else()
	set(FEELPP_ENABLE_PCH_FOR_APPLICATIONS OFF)
endif()
OPTION(FEELPP_ENABLE_RESEARCH "enable Feel++ research" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_TESTS "enable Feel++ tests" ${FEELPP_DEFAULT})
# default TimeOut for Feel++ testsuite 
SET(FEELPP_DEFAULT_TEST_TIMEOUT "300")
OPTION(FEELPP_ENABLE_DOXYGEN "enable Feel++ doxygen documentation generation" OFF)
OPTION(FEELPP_MINIMAL_CONFIGURATION "enable feel minimal configuration" OFF)
OPTION(FEELPP_MINIMAL_BUILD "Enable Feel++ minimal build" OFF)

OPTION(FEELPP_USE_GMSH_PACKAGE "Use system installed package" ON)

SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_ROOT ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Root directory for Feel++ sources")
SET(FEELPP_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Root directory for Feel++ sources")
SET(FEELPP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
SET(FEELPP_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

if ( NOT DEFINED FEELPP_INSTANTIATION_ORDER_MAX)
  # max order to instantiate
  set(FEELPP_INSTANTIATION_ORDER_MAX 3)
endif()

# if (APPLE)
#   list(REMOVE_ITEM FEELPP_LIBRARIES /usr/lib/libstdc++.dylib)
# endif()

cmake_host_system_information(RESULT NProcs QUERY NUMBER_OF_PHYSICAL_CORES)
math (EXPR NProcs2 "${NProcs}/2")
if(NProcs2 EQUAL 0)
  set(NProcs2 1)
endif()
if(NProcs2 GREATER 6)
  set(NProcs2 6)
endif()
message(STATUS "[feelpp] using ${NProcs2} cores out of ${NProcs} cores" )


SET( CMAKE_MODULE_PATH  ${FEELPP_HOME_DIR}/feelpp/cmake/modules ${FEELPP_HOME_DIR}/feelpp/cmake/machines )

include(feelpp.dirs)
include(feelpp.version)

#
# Enable testing
#
INCLUDE(CTest)
ENABLE_TESTING()

add_custom_target(check)
add_custom_target(install-testcase
  COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install-testcase
  -P cmake_install.cmake
  )

add_custom_target(check-full COMMAND "ctest")



################################################################################
# Feel
################################################################################
ADD_SUBDIRECTORY ( feelpp )

IF( FEELPP_MINIMAL_CONFIGURATION )
  set( FEELPP_ENABLE_MOR OFF )
  set( FEELPP_ENABLE_TOOLBOXES OFF )
  set( FEELPP_ENABLE_BENCHMARKS OFF )
  set( FEELPP_ENABLE_TESTS OFF )
  set( FEELPP_ENABLE_RESEARCH OFF )
  set( FEELPP_ENABLE_APPLICATIONS OFF )
  set( FEELPP_ENABLE_DOCUMENTATION ON )
  set( FEELPP_ENABLE_INSTANTIATION_MODE OFF )
  UNSET( FEELPP_INSTANTIATION_MODE CACHE )
  SET(FEELPP_MESH_MAX_ORDER "1" CACHE STRING "maximum geometrical order in templates to instantiate" FORCE )
ENDIF( FEELPP_MINIMAL_CONFIGURATION )

include(feelpp.module.toolboxes)
if ( FEELPP_HAS_TOOLBOXES )
  add_custom_target(toolboxes)
  add_subdirectory(toolboxes)
endif()

include(feelpp.module.mor)
if ( FEELPP_HAS_MOR )
  add_subdirectory(mor)
endif()

#
# Enable testing
#
IF(FEELPP_ENABLE_TESTS AND EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/testsuite)
  add_subdirectory( testsuite )
  add_dependencies(check-full testsuite)
endif()



if ( FEELPP_ENABLE_DOCUMENTATION OR FEELPP_ENABLE_DOXYGEN )
  add_dependencies(check-full doc)
  set(FEELPP_DOC_ENABLE_CAVITY_CRB off)
  ADD_SUBDIRECTORY ( doc )
endif()

if ( FEELPP_ENABLE_BENCHMARKS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks )
  ADD_SUBDIRECTORY ( benchmarks )
endif()



################################################################################
# data
################################################################################
IF(EXISTS "${CMAKE_SOURCE_DIR}/data")
  ADD_SUBDIRECTORY ( data )
ENDIF()




if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/research AND FEELPP_ENABLE_RESEARCH )
  ADD_SUBDIRECTORY ( research )
endif()

################################################################################
# generate configuration header
################################################################################
if ( CONFIG_OPTIONS )
  list(SORT CONFIG_OPTIONS)
endif()
set(FEELPP_CONFIG_OPTIONS "")
foreach(OPT ${CONFIG_OPTIONS})
  set(FEELPP_CONFIG_OPTIONS "${FEELPP_CONFIG_OPTIONS} ${OPT}")
endforeach(OPT)

################################################################################
# Installation procedure
################################################################################
include(feelpp.install)

################################################################################
# Packaging procedure
################################################################################
include(feelpp.package)

################################################################################
# Post-config message
################################################################################
include(feelpp.directive)
feelpp_print_directive()

feelpp_message( INFO ${FEELPP_MESSAGE_INFO_END} )
endif()
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
cmake_minimum_required (VERSION 2.8.7)
if ( ${CMAKE_VERSION} VERSION_GREATER 2.8.12 )
  message(STATUS "[feelpp] use new policy CMP0022" )
  cmake_policy(SET CMP0022 NEW)
endif()
if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
	message(STATUS "[feelpp] use NEW policy CMP0038" )
	cmake_policy(SET CMP0038 NEW)
	message(STATUS "[feelpp] use NEW policy CMP0042" )
	cmake_policy(SET CMP0042 NEW)
	#message(STATUS "[feelpp] use NEW policy CMP0046" )
	#cmake_policy(SET CMP0046 NEW)
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
    AND NOT cmake_build_type_tolower STREQUAL "none"
    AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are None, Debug, Release, RelWithDebInfo (case-insensitive).")
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
OPTION(FEELPP_BENCHMARK_FLAGS "enable benchmarks flags" OFF)
OPTION(FEELPP_ENABLE_DOCUMENTATION "enable Feel++ documentation (aka examples)" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_BENCHMARKS "enable Feel++ benchmarks" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_APPLICATIONS "enable Feel++ applications" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_QUICKSTART "enable Feel++ quickstart" ON)
# enable PCH only for GCC as for now
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	OPTION(FEELPP_ENABLE_PCH_FOR_APPLICATIONS "enable Pre-Compiled Headers (PCH) for Feel++ applications" OFF)
else()
	set(FEELPP_ENABLE_PCH_FOR_APPLICATIONS OFF)
endif()
OPTION(FEELPP_ENABLE_RESEARCH "enable Feel++ research" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_TESTS "enable Feel++ tests" ${FEELPP_DEFAULT})
OPTION(FEELPP_ENABLE_DOXYGEN "enable Feel++ doxygen documentation generation" OFF)
OPTION(FEELPP_MINIMAL_CONFIGURATION "enable feel minimal configuration" OFF)

OPTION(FEELPP_USE_GMSH_PACKAGE "Use system installed package" ON)

SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_ROOT ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Root directory for Feel++ sources")
SET(FEELPP_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Root directory for Feel++ sources")
SET(FEELPP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
SET(FEELPP_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")


SET( CMAKE_MODULE_PATH  ${FEELPP_HOME_DIR}/cmake/modules ${FEELPP_HOME_DIR}/cmake/machines ${FEELPP_HOME_DIR}/contrib/nt2/cmake )
add_subdirectory(cmake/modules)
add_subdirectory(cmake/machines)

# set the Feel++ version
include(feelpp.version)

# Use of standard headers for GiNaC compiler (for supermuc this option should be OFF)
option(USE_STANDARD_HEADERS_IN_GINAC_EXCOMPILER "Enable use of standard include headers in ginac-excompiler" ON)

OPTION(FEELPP_ENABLE_HOMEBREW "Enable Feel++/Homebrew support" OFF)
if (APPLE)
  set(FEELPP_ENABLE_HOMEBREW ON)
endif()

# possibly set the Feel++ configuration
include(feelpp.machines.config)

# RelWithDebInfo  is the default value for CMAKE_BUILD_TYPE
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if(APPLE)
  set(FEELPP_OS "MacOSX")
elseif(CYGWIN)
  set(FEELPP_OS "Windows")
else(APPLE)
  set(FEELPP_OS "${CMAKE_SYSTEM_NAME}")
endif(APPLE)

if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "-std=c++0x")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8.0)
     set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
     # A cosmetic change - do not work properly within emacs
     #if (CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.9.0 OR CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9.0)
     #  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always" )
     #endif()
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Default flags" FORCE)
  endif()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic  -ftemplate-depth-256 -Wno-inline" )
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Default flags" FORCE)
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fPIC" CACHE STRING "Debug flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O1" CACHE STRING "Release With Debug information flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELEASE "-g0 -O3 -DNDEBUG  -finline-functions"  CACHE STRING "Release flags" FORCE)
endif(CMAKE_COMPILER_IS_GNUCXX)

include(ProcessorCount)
ProcessorCount(NProcs)
IF(NOT CMAKE_NPROCS)
  ##if(N GREATER 1)
  math (EXPR NProcs2 '${NProcs}/2')
  if(NProcs2 EQUAL 0)
    set(NProcs2 1)
  endif()
  if(NProcs2 GREATER 6)
    set(NProcs2 6)
  endif()
else()
  set(NProcs2 "${CMAKE_NPROCS}")
endif()
message(STATUS "[feelpp] using ${NProcs2} cores out of ${NProcs} cores" )
#  set(CTEST_BUILD_FLAGS -j${N2})
#  set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N2})

# set test timeout to 300s
set(CTEST_TIMEOUT "20")

#else()
#  OPTION(FEELPP_ENABLE_MPI_MODE "Instantiation mode" OFF )
#endif()


if ( FEELPP_BENCHMARK_FLAGS )
  set(CMAKE_BUILD_TYPE Release )
  set(GCC_PARAM_INLINE_UNIT_GROWTH 150)
  set(GCC_PARAM_MAX_INLINE_INSNS_SINGLE 500)
  set(GCC_PARAM_LARGE_FUNCTION_GROWTH 600)
  #  SET(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wshadow -Woverloaded-virtual -std=c++0x -O3 -DNDEBUG --param max-inline-recursive-depth=256 --param max-gcse-memory=8000 --param max-inline-insns-single=${GCC_PARAM} --param inline-unit-growth=${GCC_PARAM} --param large-unit-insns=${GCC_PARAM} --param large-function-growth=${GCC_PARAM} --param large-function-insns=${GCC_PARAM} " CACHE STRING "Benchmarks Release flags" FORCE)
  #  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -DNDEBUG --param max-inline-recursive-depth=256 --param max-gcse-memory=8000 --param max-inline-insns-single=${GCC_PARAM} --param inline-unit-growth=${GCC_PARAM} --param large-unit-insns=${GCC_PARAM} --param large-function-growth=${GCC_PARAM} --param large-function-insns=${GCC_PARAM} " CACHE STRING "Benchmarks Release flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -DNDEBUG --param max-inline-insns-single=${GCC_PARAM_MAX_INLINE_INSNS_SINGLE} --param inline-unit-growth=${GCC_PARAM_INLINE_UNIT_GROWTH} --param large-function-growth=${GCC_PARAM_LARGE_FUNCTION_GROWTH} " CACHE STRING "Benchmarks Release flags" FORCE)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  option(FEELPP_ENABLE_SSE2 "Enable/Disable SSE2 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    message(STATUS "[feelpp] Enabling SSE2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE3 "Enable/Disable SSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    message(STATUS "[feelpp] Enabling SSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSSE3 "Enable/Disable SSSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
    message(STATUS "[feelpp] Enabling SSSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_1 "Enable/Disable SSE4.1 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    message(STATUS "[feelpp] Enabling SSE4.1 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_2 "Enable/Disable SSE4.2 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    message(STATUS "[feelpp] Enabling SSE4.2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_ALTIVEC "Enable/Disable AltiVec in tests/examples" OFF)
  if(FEELPP_ENABLE_ALTIVEC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maltivec -mabi=altivec")
    message(STATUS "[feelpp] Enabling AltiVec in tests/examples")
  endif()

  option(FEELPP_ENABLE_NEON "Enable/Disable Neon in tests/examples" OFF)
  if(FEELPP_ENABLE_NEON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=softfp -mfpu=neon -mcpu=cortex-a8")
    message(STATUS "[feelpp] Enabling NEON in tests/examples")
  endif()
endif(CMAKE_COMPILER_IS_GNUCXX)
#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(FeelGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)

include(feelpp.extra.warnings)
include(feelpp.extra.astyle)
add_definitions(${FEELPP_FLAGS})

if("${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}" MATCHES "[^ ]")
  message(STATUS "[feelpp] [feel++] Global flags: ${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}")
endif()

message(STATUS "[feelpp] [feel++] Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "[feelpp] [feel++] Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")

INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)
INCLUDE(ParseArguments)

# MACRO(CAR var)
#   SET(${var} ${ARGV1})
# ENDMACRO(CAR)

# MACRO(CDR var junk)
#   SET(${var} ${ARGN})
# ENDMACRO(CDR)
INCLUDE(CheckTypeSize)
CHECK_TYPE_SIZE(int SIZE_INT )
CHECK_TYPE_SIZE(uint SIZE_UINT )
CHECK_TYPE_SIZE(size_t SIZE_SIZE_T )
CHECK_TYPE_SIZE(long SIZE_LONG )
CHECK_TYPE_SIZE(float SIZE_FLOAT )
CHECK_TYPE_SIZE(double SIZE_DOUBLE )
CHECK_TYPE_SIZE("long double" SIZE_LONG_DOUBLE)
MESSAGE(STATUS "[feelpp] SIZE_INT=${SIZE_INT}")
MESSAGE(STATUS "[feelpp] SIZE_UINT=${SIZE_UINT}")
MESSAGE(STATUS "[feelpp] SIZE_SIZE_T=${SIZE_SIZE_T}")
MESSAGE(STATUS "[feelpp] SIZE_LONG=${SIZE_LONG}")
MESSAGE(STATUS "[feelpp] SIZE_FLOAT=${SIZE_FLOAT}")
MESSAGE(STATUS "[feelpp] SIZE_DOUBLE=${SIZE_DOUBLE}")
MESSAGE(STATUS "[feelpp] SIZE_LONG_DOUBLE=${SIZE_LONG_DOUBLE}")

# Find Feel++ or at least it dependencies
FIND_PACKAGE(Feel++ REQUIRED)

include( feelpp.macros )

# Load documented project for Doxygen
# it assumes that a project is valid if :
#   there is README.md file in the directory
#   that the directory path include "research"
# for each detected project it creates an FEELPP_ENABLE_${PROJECT} and FEELPP_ENABLE_${PROJECT}_DOCUMENTATION
# option which are set to OFF by default
MESSAGE(STATUS "[feelpp] Feel++ projects :")
feelpp_list_subdir(projects ${CMAKE_CURRENT_SOURCE_DIR}/research)
foreach(project ${projects})
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/research/${project}/README.md )
    string(TOUPPER ${project} PROJECT)
    option(FEELPP_ENABLE_RESEARCH_${PROJECT} "Enable/Disable research/${project} " OFF)
    option(FEELPP_ENABLE_${PROJECT}_DOCUMENTATION "Enable/Disable documentation in research/${project}" OFF)
		SET(FEELPP_ENABLED_PROJECTS "${FEELPP_ENABLED_PROJECTS} ${project}" )
		MESSAGE(STATUS "[feelpp] ${project}\t-- Compile : ${FEELPP_ENABLE_RESEARCH_${PROJECT}} ")
		MESSAGE(STATUS "        \t\t-- Doc :     ${FEELPP_ENABLE_${PROJECT}_DOCUMENTATION}" )
  endif()
endforeach()

################################################################################
# data
################################################################################
IF(EXISTS "${CMAKE_SOURCE_DIR}/data")
  ADD_SUBDIRECTORY ( data )
ENDIF()

if (APPLE)
  list(REMOVE_ITEM FEELPP_LIBRARIES /usr/lib/libstdc++.dylib)
endif()



################################################################################
# Feel
################################################################################
if(NOT FEELPP_ENABLE_BUILD_STATIC)
	set(CMAKE_C++_CREATE_SHARED_LIBRARY 1)
endif()
ADD_SUBDIRECTORY ( feel )

IF( FEELPP_MINIMAL_CONFIGURATION )
  set( FEELPP_ENABLE_BENCHMARKS OFF )
  set( FEELPP_ENABLE_TESTS OFF )
  set( FEELPP_ENABLE_RESEARCH OFF )
  set( FEELPP_ENABLE_APPLICATIONS OFF )
  set( FEELPP_ENABLE_DOCUMENTATION ON )
  set( FEELPP_ENABLE_INSTANTIATION_MODE OFF )
  UNSET( FEELPP_INSTANTIATION_MODE CACHE )
  SET(FEELPP_MESH_MAX_ORDER "1" CACHE STRING "maximum geometrical order in templates to instantiate" FORCE )
ENDIF( FEELPP_MINIMAL_CONFIGURATION )

option(FEELPP_ENABLE_APPLICATIONS_CRB "Enable CRB applications in Feel++" ON)
if (ANN_FOUND AND GLPK_FOUND)
  #set( FEELPP_ENABLE_APPLICATIONS_CRB ON )
  #set( FEELPP_ENABLE_OPENTURNS ON)
else()
  set( FEELPP_ENABLE_APPLICATIONS_CRB OFF )
endif()

#
# Enable testing
#
INCLUDE(CTest)
ENABLE_TESTING()

add_custom_target(check)

add_custom_target(check-full COMMAND "ctest")
IF(FEELPP_ENABLE_TESTS )
  add_dependencies(check-full testsuite)
endif()
if ( FEELPP_ENABLE_DOCUMENTATION OR FEELPP_ENABLE_DOXYGEN )
  add_dependencies(check-full doc)
  set(FEELPP_DOC_ENABLE_CAVITY_CRB off)
  ADD_SUBDIRECTORY ( doc )
endif()
#add_dependencies(check benchmarks)
#add_dependencies(check examples)
IF(FEELPP_ENABLE_TESTS AND EXISTS  ${CMAKE_CURRENT_SOURCE_DIR}/testsuite)
  add_subdirectory( testsuite )
endif()

if ( FEELPP_ENABLE_BENCHMARKS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks )
  ADD_SUBDIRECTORY ( benchmarks )
endif()

if ( FEELPP_ENABLE_APPLICATIONS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/applications)
  ADD_SUBDIRECTORY ( applications )
elseif (EXISTS ${FEELPP_SOURCE_DIR}/applications/opus )
  ADD_SUBDIRECTORY ( applications/opus )
endif()

if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/research AND FEELPP_ENABLE_RESEARCH )
  ADD_SUBDIRECTORY ( research )
endif()

include(feelpp.module.quickstart)


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

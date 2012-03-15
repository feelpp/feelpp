# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
cmake_minimum_required (VERSION 2.6)

# guard against in-source builds

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# guard against bad build-type strings

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release"
   AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()


#SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
#SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

project (Feel++ C CXX Fortran)

macro(set_config_option VARNAME STRING)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "Found " ${STRING})
endmacro(set_config_option)


OPTION(FEELPP_ENABLE_SVN "enable Feel++ looking up for svn information" OFF)
SET(FEELPP_SCM "svn")
FIND_PACKAGE(Subversion)
IF(FEELPP_ENABLE_SVN AND ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn ) )
#IF ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Current revision is ${Project_WC_REVISION}")
  Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Last changed log is ${Project_LAST_CHANGED_LOG}")
ELSE()
  FIND_PACKAGE(Git)
  if(FEELPP_ENABLE_SVN AND GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git )

    execute_process(
      COMMAND git rev-parse --verify -q --short git-svn
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_COMMIT
      ERROR_VARIABLE GIT_COMMIT_error
      RESULT_VARIABLE GIT_COMMIT_result)
    if ( NOT ${GIT_COMMIT_result} EQUAL 0 )
      MESSAGE(SEND_ERROR "Command \"git rev-parse --verify -q --short git-svn\" failed with output:\n${GIT_COMMIT_error}")
    else()
      STRING(STRIP "${GIT_COMMIT}" GIT_COMMIT )
      MESSAGE(STATUS "Git commit: ${GIT_COMMIT}")
    endif()
    execute_process(
      COMMAND git show --raw ${GIT_COMMIT}
#      COMMAND "grep git-svn-id"
#      COMMAND sed -re "s/^\\s*git-svn-id: .*@([0-9]+).*$/\\1/"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE SVN_COMMIT
      ERROR_VARIABLE Project_WC_REVISION_error
      RESULT_VARIABLE Project_WC_REVISION_result)
    if ( NOT ${Project_WC_REVISION_result} EQUAL 0 )
      MESSAGE(SEND_ERROR "Command \"git show --raw ${GIT_COMMIT}\" failed with output:\n${Project_WC_REVISION_error}")
    endif()
    string(REGEX REPLACE ".*git-svn-id: .*trunk@([0-9]+) .*" "\\1"Project_WC_REVISION ${SVN_COMMIT})
    MESSAGE(STATUS "Current SVN(git) revision is ${Project_WC_REVISION}")
  endif()
ENDIF()

set(FEELPP_VERSION_MAJOR "0")
set(FEELPP_VERSION_MINOR "91")
set(FEELPP_VERSION_MICRO "0")
set(FEELPP_REVISION "0" )
set(FEELPP_BUILDID "0" )
if ((FEELPP_ENABLE_SVN AND Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) OR
    (FEELPP_ENABLE_SVN AND GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git ))
  set(FEELPP_VERSION_EXTRA "~${FEELPP_SCM}${Project_WC_REVISION}")
  set(FEELPP_REVISION "${Project_WC_REVISION}")
  set(FEELPP_BUILDID "${Project_WC_REVISION}")
  file(WRITE "${CMAKE_SOURCE_DIR}/SVNREVISION" "${Project_WC_REVISION}")
else()
  set(FEELPP_VERSION_EXTRA "")
  if ( EXISTS "${CMAKE_SOURCE_DIR}/SVNREVISION" )
    file(READ "${CMAKE_SOURCE_DIR}/SVNREVISION" FEELPP_REVISION)
    file(READ "${CMAKE_SOURCE_DIR}/SVNREVISION" FEELPP_BUILDID)
  endif()
endif()
set(FEELPP_VERSION "(((${FEELPP_VERSION_MAJOR}) << 16) | ((${FEELPP_VERSION_MINOR}) << 9) | (${FEELPP_VERSION_MICRO}))" )
set(FEELPP_VERSION_STRING "${FEELPP_VERSION_MAJOR}.${FEELPP_VERSION_MINOR}.${FEELPP_VERSION_MICRO}${FEELPP_VERSION_EXTRA}")
set(FEELPP_SHARED_VERSION "1.0.0" )
set(FEELPP_SHARED_SOVERSION "1" )

SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

SET( CMAKE_MODULE_PATH  ${FEELPP_HOME_DIR}/cmake/modules )
add_subdirectory(cmake/modules)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
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
  SET(CMAKE_CXX_FLAGS "-std=c++0x -pedantic  -ftemplate-depth-256 -Wno-inline" CACHE STRING "Default flags" FORCE)
  SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x  -pedantic -g  -ftemplate-depth-256 -Wno-inline" CACHE STRING "Debug flags" FORCE)
  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x  -pedantic -g0 -O2 -DNDEBUG  -ftemplate-depth-256 -finline-functions -Wno-inline"  CACHE STRING "Release flags" FORCE)
endif(CMAKE_COMPILER_IS_GNUCXX)

OPTION(FEELPP_ENABLE_MOVE_SEMANTICS "enable move semantics(elision)" ON )
MARK_AS_ADVANCED(FEELPP_ENABLE_MOVE_SEMANTICS)
IF ( FEELPP_ENABLE_MOVE_SEMANTICS )
  SET( BOOST_UBLAS_MOVE_SEMANTICS 1 CACHE STRING "Enable Boost Ublas move semantics" FORCE )
  ADD_DEFINITIONS( -DBOOST_UBLAS_MOVE_SEMANTICS )
ENDIF( FEELPP_ENABLE_MOVE_SEMANTICS )


OPTION(ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
MARK_AS_ADVANCED(ENABLE_INSTANTIATION_MODE)
IF ( ENABLE_INSTANTIATION_MODE )
  SET( FEELPP_INSTANTIATION_MODE 1 )
ENDIF()

OPTION(FEELPP_ENABLE_MPI_MODE "Instantiation mode" OFF )
IF ( FEELPP_ENABLE_MPI_MODE )
  SET( FEELPP_ENABLE_MPI_MODE 1 )
ENDIF()

OPTION(FEELPP_BENCHMARK_FLAGS "enable benchmarks flags" OFF)
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
    message(STATUS "Enabling SSE2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE3 "Enable/Disable SSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    message(STATUS "Enabling SSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSSE3 "Enable/Disable SSSE3 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSSE3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
    message(STATUS "Enabling SSSE3 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_1 "Enable/Disable SSE4.1 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    message(STATUS "Enabling SSE4.1 in tests/examples")
  endif()

  option(FEELPP_ENABLE_SSE4_2 "Enable/Disable SSE4.2 in tests/examples" OFF)
  if(FEELPP_ENABLE_SSE4_2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
    message(STATUS "Enabling SSE4.2 in tests/examples")
  endif()

  option(FEELPP_ENABLE_ALTIVEC "Enable/Disable AltiVec in tests/examples" OFF)
  if(FEELPP_ENABLE_ALTIVEC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -maltivec -mabi=altivec")
    message(STATUS "Enabling AltiVec in tests/examples")
  endif()

  option(FEELPP_ENABLE_NEON "Enable/Disable Neon in tests/examples" OFF)
  if(FEELPP_ENABLE_NEON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=softfp -mfpu=neon -mcpu=cortex-a8")
    message(STATUS "Enabling NEON in tests/examples")
  endif()
endif(CMAKE_COMPILER_IS_GNUCXX)
#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(FeelGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)

include(feelpp.extra.warnings)
add_definitions(${FEELPP_FLAGS})

if("${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}" MATCHES "[^ ]")
  message(STATUS "[feel++] Global flags: ${CMAKE_CXX_FLAGS} ${FEELPP_FLAGS}")
endif()

message(STATUS "[feel++] Debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "[feel++] Release flags: ${CMAKE_CXX_FLAGS_RELEASE}")

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
MESSAGE(STATUS "SIZE_INT=${SIZE_INT}")
MESSAGE(STATUS "SIZE_UINT=${SIZE_UINT}")
MESSAGE(STATUS "SIZE_SIZE_T=${SIZE_SIZE_T}")
MESSAGE(STATUS "SIZE_LONG=${SIZE_LONG}")
MESSAGE(STATUS "SIZE_FLOAT=${SIZE_FLOAT}")
MESSAGE(STATUS "SIZE_DOUBLE=${SIZE_DOUBLE}")
MESSAGE(STATUS "SIZE_LONG_DOUBLE=${SIZE_LONG_DOUBLE}")


# Find Feel++ or at least it dependencies
FIND_PACKAGE(Feel++ REQUIRED)

option(FEELPP_ENABLE_BAMG "Enable Bamg mesh generator in Feel++" ON)
if(FEELPP_ENABLE_BAMG)
  add_subdirectory(contrib/bamg)
  include_directories(contrib/bamg contrib/bamg/bamglib)
  set_config_option(FEELPP_HAVE_BAMG "Bamg")
endif(FEELPP_ENABLE_BAMG)


#
# Python
#
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)

#
# data
#
IF(NOT APPLE)
ADD_SUBDIRECTORY ( data )
ENDIF(NOT APPLE)

#
# Feel
#
ADD_SUBDIRECTORY ( feel )




SET(FEELPP_MESH_MAX_ORDER "5" CACHE STRING "maximum geometrical order in templates to instantiate" )

OPTION(FEELPP_ENABLE_DOCUMENTATION "enable Feel++ documentation" ON)
OPTION(FEELPP_ENABLE_BENCHMARKS "enable Feel++ benchmarks" OFF)
OPTION(FEELPP_ENABLE_EXAMPLES "enable Feel++ examples" ON)
OPTION(FEELPP_ENABLE_APPLICATIONS "enable Feel++ applications" ON)
OPTION(FEELPP_ENABLE_RESEARCH "enable Feel++ research" ON)
OPTION(FEELPP_ENABLE_TESTS "enable Feel++ tests" ON)

OPTION(FEELPP_MINIMAL_CONFIGURATION "enable feel minimal configuration" OFF)
IF( FEELPP_MINIMAL_CONFIGURATION )
  set( BUILD_TESTING OFF )
  set( FEELPP_ENABLE_EXAMPLES OFF )
  set( FEELPP_ENABLE_BENCHMARKS OFF )
  set( FEELPP_ENABLE_TESTS OFF )
  set( FEELPP_ENABLE_RESEARCH OFF )
  set( FEELPP_ENABLE_APPLICATIONS OFF )
  set( FEELPP_ENABLE_DOCUMENTATION ON )
  set( ENABLE_INSTANTIATION_MODE OFF )
  UNSET( FEELPP_INSTANTIATION_MODE CACHE )
  SET(FEELPP_MESH_MAX_ORDER "1" CACHE STRING "maximum geometrical order in templates to instantiate" FORCE )
ENDIF( FEELPP_MINIMAL_CONFIGURATION )

option(FEELPP_ENABLE_APPLICATIONS_CRB "Enable CRB applications in Feel++" OFF)
if (ANN_FOUND AND GLPK_FOUND)
  set( FEELPP_ENABLE_APPLICATIONS_CRB ON )
  set( FEELPP_ENABLE_OPENTURNS ON)
endif()

#
# Enable testing
#
INCLUDE(CTest)
ENABLE_TESTING()
add_custom_target(check COMMAND "ctest")
add_dependencies(check testsuite)
add_dependencies(check doc)
#add_dependencies(check benchmarks)
#add_dependencies(check examples)

if ( FEELPP_ENABLE_DOCUMENTATION )
  ADD_SUBDIRECTORY ( doc )
endif ( FEELPP_ENABLE_DOCUMENTATION )

if ( FEELPP_ENABLE_BENCHMARKS )
  ADD_SUBDIRECTORY ( benchmarks )
endif()

if ( FEELPP_ENABLE_EXAMPLES )
 ADD_SUBDIRECTORY ( examples )
endif()

if ( FEELPP_ENABLE_APPLICATIONS )
 ADD_SUBDIRECTORY ( applications )
elseif (EXISTS ${FEELPP_SOURCE_DIR}/applications/opus )
  ADD_SUBDIRECTORY ( applications/opus )
endif()


if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/research AND FEELPP_ENABLE_RESEARCH )
  ADD_SUBDIRECTORY ( research )
endif()




IF(FEELPP_ENABLE_TESTS )

  macro(feel_add_test testname)
    set(targetname test_${testname})
    set(filename test_${testname}.cpp)
    add_executable(${targetname} ${filename})
    target_link_libraries(${targetname} ${FEELPP_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
    set_property(TARGET ${targetname} PROPERTY LABELS testsuite)
    add_dependencies(testsuite ${targetname})
    if ( APPLE )
      add_test(
        NAME test_${testname}
        COMMAND ${targetname}  # --log_level=message seems to trigger a segfault
        )
    else( APPLE )
      add_test(
        NAME test_${testname}
        COMMAND ${targetname} --log_level=message
        )
    endif ( APPLE )
    set_property(TEST ${targetname} PROPERTY LABELS testsuite)
  endmacro(feel_add_test)

  add_subdirectory( testsuite )

endif()

# generate configuration header
if ( CONFIG_OPTIONS )
  list(SORT CONFIG_OPTIONS)
endif()
set(FEELPP_CONFIG_OPTIONS "")
foreach(OPT ${CONFIG_OPTIONS})
  set(FEELPP_CONFIG_OPTIONS "${FEELPP_CONFIG_OPTIONS} ${OPT}")
endforeach(OPT)

set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEELPP_DATADIR ${CMAKE_INSTALL_PREFIX}/share/feel )
CONFIGURE_FILE(feelconfig.h.in feel/feelconfig.h  @ONLY)
CONFIGURE_FILE(feelinfo.h.in feel/feelinfo.h  @ONLY)

#
# Packaging
#
INCLUDE(InstallRequiredSystemLibraries)


foreach(includedir feelcore feelalg feelmesh feelpoly feelfilters feeldiscr feelvf feelmaterial feelsystem )
  FILE(GLOB files "feel/${includedir}/*.hpp")
  FILE(GLOB cppfiles "feel/${includedir}/*.cpp")
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/feel/${includedir} COMPONENT Devel)
endforeach()
FILE(GLOB files "feel/*.hpp")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
# gmm
FILE(GLOB files "contrib/gmm/include/*.h")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
# feel++ config headers
FILE(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/feel/*.h")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)

message(STATUS "apps: ${FEELPP_INSTALL_APPS}")
INSTALL(FILES ${FEELPP_INSTALL_APPS} DESTINATION bin/ COMPONENT Bin)

# documentation and examples
if ( FEELPP_ENABLE_DOCUMENTATION )
#  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/feel_get_tutorial.sh DESTINATION bin COMPONENT Doc)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/manual/feel-manual.pdf DESTINATION share/doc/feel COMPONENT Doc)
  if ( NOT APPLE )
    if ( EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc/api/html )
      install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/api/html DESTINATION share/doc/feel COMPONENT Doc
        PATTERN ".svn" EXCLUDE PATTERN ".git" EXCLUDE )
    elseif( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/api/html )
      install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/api/html DESTINATION share/doc/feel COMPONENT Doc
        PATTERN ".svn" EXCLUDE PATTERN ".git" EXCLUDE )
    endif()
  endif(NOT APPLE)
  FILE(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/doc/manual/*.*pp")
  INSTALL(FILES ${examples} DESTINATION share/doc/feel/examples/ COMPONENT Doc)
endif()


##
## Archive generation using cpack
##
if (UNIX)
  execute_process(
    COMMAND uname -m
    OUTPUT_VARIABLE FEELPP_SYSTEM_MACHINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE FEELPP_SYSTEM_MACHINE_error
    RESULT_VARIABLE FEELPP_SYSTEM_MACHINE_result)
endif()
SET(CPACK_PACKAGE_NAME "feel++")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FEELPP_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FEELPP_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FEELPP_VERSION_MICRO}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "feel")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
SET(CPACK_SYSTEM_NAME "${FEELPP_OS}-${FEELPP_SYSTEM_MACHINE}")

OPTION(FEELPP_ENABLE_CPACK_OPUS "Enable OPUS packaging (if available) in CPack along with Feel++" ON )
SET(CPACK_PACKAGE_NAME "feel++")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++ (Finite Element method Embedded Library and language in C++)")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${FEELPP_VERSION_MAJOR}.${FEELPP_VERSION_MINOR}.${FEELPP_VERSION_MICRO}${FEELPP_VERSION_EXTRA}")


SET(CPACK_SOURCE_STRIP_FILES "")
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git/;\\\\.gitignore;/\\\\.svn;"
  "/admin/;/tools/;/Templates/;"
  "/auto/;/auto.*/;"
  "/TAGS;/#.*;/.*~;/.cvsignore;/.bzrignore;/work/"
  "${PROJECT_SOURCE_DIR}/benchmarks/turek/"
  "${PROJECT_SOURCE_DIR}/benchmarks/stokes/"
  "${PROJECT_SOURCE_DIR}/benchmarks/ethiersteinman/"
  "${PROJECT_SOURCE_DIR}/benchmarks/kovasznay/"
  "${PROJECT_SOURCE_DIR}/doc/poster/"
  "${PROJECT_SOURCE_DIR}/doc/manual/pdfs/"
  "${PROJECT_SOURCE_DIR}/doc/figures/backgrounds/"
  "${PROJECT_SOURCE_DIR}/doc/figures/logos/"
  "${PROJECT_SOURCE_DIR}/examples/fluid/"
  "${PROJECT_SOURCE_DIR}/examples/levelset/"
  "${PROJECT_SOURCE_DIR}/examples/pbeq/"
  "${PROJECT_SOURCE_DIR}/research/"
  "/*.tar.gz;*.tar.bz2;*.deb;obj-x86_64-linux-gnu/;/opt/"
  "/TAGS;/#.*;/.*~;/.cvsignore;/.bzrignore"
  "*.aux;*.log;*.bbl;*.idx;*.ist;*.out;*.blg;OpusManualBenchmarkEADSUJF.pdf"
  "${PROJECT_SOURCE_DIR}/applications/opus.old"
  "${PROJECT_SOURCE_DIR}/applications/opus/"
   "${PROJECT_SOURCE_DIR}/applications/opus/debian/opus/"
  "${PROJECT_SOURCE_DIR}/applications/opus/debian/source/"
  "${PROJECT_SOURCE_DIR}/applications/opus/doc/"
  "${PROJECT_SOURCE_DIR}/applications/opus/scripts"
  )

#if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/applications/opus/ )
if ( EXISTS ${FEELPP_SOURCE_DIR}/applications/opus/ AND NOT FEELPP_ENABLE_CPACK_OPUS )
#  set(CPACK_SOURCE_IGNORE_FILES "${PROJECT_SOURCE_DIR}/applications/opus/;${CPACK_SOURCE_IGNORE_FILES}" )
endif()

include( CPack )


MESSAGE(STATUS "================================================================================")
MESSAGE(STATUS "FEELPP_VERSION_MAJOR : ${FEELPP_VERSION_MAJOR}")
MESSAGE(STATUS "FEELPP_VERSION_MINOR : ${FEELPP_VERSION_MINOR}")
MESSAGE(STATUS "FEELPP_VERSION_MICRO : ${FEELPP_VERSION_MICRO}")
MESSAGE(STATUS "FEELPP_REVISON : ${FEELPP_REVISION}")
MESSAGE(STATUS "FEELPP_BUILDID : ${FEELPP_BUILDID}")
MESSAGE(STATUS "")
MESSAGE(STATUS "Feel++ Modules :")
MESSAGE(STATUS "  Documentation: ${FEELPP_ENABLE_DOCUMENTATION}")
MESSAGE(STATUS "      Testsuite: ${FEELPP_ENABLE_TESTS}")
MESSAGE(STATUS "   Applications: ${FEELPP_ENABLE_APPLICATIONS}")
MESSAGE(STATUS "       Examples: ${FEELPP_ENABLE_EXAMPLES}")
MESSAGE(STATUS "     Benchmarks: ${FEELPP_ENABLE_BENCHMARKS}")
MESSAGE(STATUS "       Research: ${FEELPP_ENABLE_RESEARCH}")
MESSAGE(STATUS "")
MESSAGE(STATUS "Feel++ Configuration:")
MESSAGE(STATUS "         CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "       CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "          CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "  CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
MESSAGE(STATUS "    CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE(STATUS "")
MESSAGE(STATUS "          FEELPP_ENABLE_SVN: ${FEELPP_ENABLE_SVN}")
MESSAGE(STATUS "      FEELPP_MESH_MAX_ORDER: ${FEELPP_MESH_MAX_ORDER}")
MESSAGE(STATUS "  FEELPP_INSTANTIATION_MODE: ${FEELPP_INSTANTIATION_MODE}")
MESSAGE(STATUS "     FEELPP_ENABLED_OPTIONS: ${FEELPP_ENABLED_OPTIONS}")
MESSAGE(STATUS "================================================================================")

string(TOLOWER "${CMAKE_GENERATOR}" cmake_generator_tolower)
if(cmake_generator_tolower MATCHES "makefile")
  message(STATUS "Some things you can do now with Feel++:")
  MESSAGE(STATUS "================================================================================")
  message(STATUS "Command        |   Description")
  MESSAGE(STATUS "===============|================================================================")
  message(STATUS "make install   | Install to ${CMAKE_INSTALL_PREFIX}. To change that:")
  message(STATUS "               |     cmake . -DCMAKE_INSTALL_PREFIX=yourpath")
  message(STATUS "               |   Feel++ headers will then be installed to:")
  message(STATUS "               |     ${INCLUDE_INSTALL_DIR}")
  message(STATUS "make doc       | Generate the manual applications, manual and the API ")
  message(STATUS "               | documentation, requires Doxygen & LaTeX")
  message(STATUS "make examples  | Generate the examples ")
  message(STATUS "make benchmarks| Generate the benchmarks(Warning: takes a long time!) ")
  message(STATUS "make check     | Build and run the unit-tests.")
  MESSAGE(STATUS "================================================================================")
endif()

message(STATUS "")

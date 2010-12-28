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

SET(CMAKE_CXX_FLAGS "-std=c++0x  -mssse3" CACHE STRING "Default flags" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x -g  -mssse3" CACHE STRING "Debug flags" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O2  -mssse3 -DNDEBUG" CACHE STRING "Release flags" FORCE)


LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG)
LIST(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELEASE)

#SET( CMAKE_CXX_FLAGS "-pipe -Wall -O2 ")
#SET( CMAKE_C_FLAGS "-pipe -Wall -O2")

project (Feel C CXX Fortran)

OPTION(FEEL_ENABLE_SVN "enable Feel++ looking up for svn information" ON)
SET(FEEL_SCM "svn")
FIND_PACKAGE(Subversion)
#IF(FEEL_ENABLE_SVN AND ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn ) )
IF ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Current revision is ${Project_WC_REVISION}")
  Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Last changed log is ${Project_LAST_CHANGED_LOG}")
ELSE()
  FIND_PACKAGE(Git)
  if(GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git )

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

set(FEEL_VERSION_MAJOR "0")
set(FEEL_VERSION_MINOR "90")
set(FEEL_VERSION_MICRO "0")
if ((Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) OR
    (GIT_FOUND AND  EXISTS ${PROJECT_SOURCE_DIR}/.git ))
  set(FEEL_VERSION_EXTRA "~${FEEL_SCM}${Project_WC_REVISION}")
  set(FEEL_REVISION "${Project_WC_REVISION}")
  set(FEEL_BUILDID "${Project_WC_REVISION}")
  file(WRITE "${CMAKE_SOURCE_DIR}/SVNREVISION" "${Project_WC_REVISION}")
else()
  set(FEEL_VERSION_EXTRA "")
  if ( EXISTS "${CMAKE_SOURCE_DIR}/SVNREVISION" )
    file(READ "${CMAKE_SOURCE_DIR}/SVNREVISION" FEEL_REVISION)
    file(READ "${CMAKE_SOURCE_DIR}/SVNREVISION" FEEL_BUILDID)
  else()
    set(FEEL_REVISION "0")
    set(FEEL_BUILDID "0")
  endif()
endif()
set(FEEL_VERSION "(((${FEEL_VERSION_MAJOR}) << 16) | ((${FEEL_VERSION_MINOR}) << 9) | (${FEEL_VERSION_MICRO}))" )
set(FEEL_VERSION_STRING "${FEEL_VERSION_MAJOR}.${FEEL_VERSION_MINOR}.${FEEL_VERSION_MICRO}${FEEL_VERSION_EXTRA}")
set(FEEL_SHARED_VERSION "1.0.0" )
set(FEEL_SHARED_SOVERSION "1" )

SET(FEEL_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEEL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEEL_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

SET( CMAKE_MODULE_PATH  ${FEEL_HOME_DIR}/cmake/modules )

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

OPTION(ENABLE_MOVE_SEMANTICS "enable move semantics(elision)" ON )
MARK_AS_ADVANCED(ENABLE_MOVE_SEMANTICS)
IF ( ENABLE_MOVE_SEMANTICS )
  SET( BOOST_UBLAS_MOVE_SEMANTICS 1 CACHE STRING "Enable Boost Ublas move semantics" FORCE )
  ADD_DEFINITIONS( -DBOOST_UBLAS_MOVE_SEMANTICS )
ENDIF( ENABLE_MOVE_SEMANTICS)

OPTION(ENABLE_INSTANTIATION_MODE "Instantiation mode" ON )
MARK_AS_ADVANCED(ENABLE_INSTANTIATION_MODE)
IF ( ENABLE_INSTANTIATION_MODE )
  SET( FEEL_INSTANTIATION_MODE 1 )
ENDIF()

OPTION(FEEL_BENCHMARK_FLAGS "enable benchmarks flags" OFF)
if ( FEEL_BENCHMARK_FLAGS )
  set(CMAKE_BUILD_TYPE Release )
  SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++0x -O3 -mssse3 -DNDEBUG --param max-inline-recursive-depth=32 --param max-gcse-passes=3 --param max-gcse-memory=8000 --param max-inline-insns-single=2000" CACHE STRING "Release flags" FORCE)
endif()


#INCLUDE(PackageArchGlobalMacros)
# INCLUDE(FeelGlobalMacros)
#INCLUDE(AdvancedSet)
#INCLUDE(AdvancedOption)


INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckLibraryExists)
# INCLUDE(parse_arguments)

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


OPTION(FEEL_ENABLE_TBB "enable feel++ TBB support" OFF)
if ( FEEL_ENABLE_TBB )
  FIND_PACKAGE(TBB)
  IF ( TBB_FOUND )
    INCLUDE_DIRECTORIES( ${TBB_INCLUDE_DIR} )
    SET(FEEL_LIBRARIES ${TBB_LIBRARIES})
  ENDIF (TBB_FOUND )
endif()

FIND_PACKAGE(OpenMP)

# on APPLE enfore the use of macports openmpi version
if ( APPLE )
set(MPI_COMPILER /opt/local/bin/mpic++)
set(MPI_LIBRARY "MPI_LIBRARY-NOTFOUND" )
endif( APPLE )
FIND_PACKAGE(MPI)
IF ( MPI_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  SET( HAVE_MPI 1 )
  SET( HAVE_MPI_H 1 )
  ADD_DEFINITIONS( -DHAVE_MPI=1 -DHAVE_MPI_H=1 )
  SET(FEEL_BOOST_MPI mpi)
  SET(FEEL_LIBRARIES ${MPI_LIBRARIES} ${FEEL_LIBRARIES})
ENDIF()



Check_Include_File_CXX(dlfcn.h HAVE_DLFCN_H)
if ( HAVE_DLFCN_H )
  add_definitions(-DHAVE_DLFCN_H)
endif()
CHECK_LIBRARY_EXISTS (dl dlopen "" HAVE_LIBDL)
IF (HAVE_LIBDL)
  set(DL_LIBRARY dl)
  SET (HAVE_DLOPEN 1)
  add_definitions(-DHAVE_DLOPEN)
  SET(FEEL_LIBRARIES ${DL_LIBRARY} ${FEEL_LIBRARIES})
ELSE ()
  CHECK_FUNCTION_EXISTS (dlopen HAVE_DLOPEN)
ENDIF (HAVE_LIBDL)

#
# Xml2
#
FIND_PACKAGE(LibXml2 REQUIRED)
SET(FEEL_LIBRARIES  ${LIBXML2_LIBRARIES} ${FEEL_LIBRARIES})
#
# Blas and Lapack
#
if (APPLE)
        FIND_PACKAGE(LAPACK )
else (APPLE)
        FIND_PACKAGE(LAPACK REQUIRED)
endif (APPLE)
SET(FEEL_LIBRARIES  ${LAPACK_LIBRARIES} ${FEEL_LIBRARIES})

FIND_PACKAGE(Boost COMPONENTS date_time filesystem system program_options unit_test_framework signals  ${FEEL_BOOST_MPI} regex  serialization)
set(Boost_ADDITIONAL_VERSIONS "1.39" "1.40" "1.41" "1.42" "1.43" "1.44" "1.45" )
set( BOOST_PARAMETER_MAX_ARITY 15 )
add_definitions( -DBOOST_PARAMETER_MAX_ARITY=${BOOST_PARAMETER_MAX_ARITY} -DBOOST_TEST_DYN_LINK )
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

SET(FEEL_LIBRARIES ${Boost_LIBRARIES} ${FEEL_LIBRARIES})

#
# MadLib
#
FIND_PACKAGE(MadLib)
if ( MADLIB_FOUND )
  INCLUDE_DIRECTORIES(${MadLib_INCLUDE_DIR})
  SET (HAVE_MADLIB_H 1)
  LINK_DIRECTORIES(${MadLib_LIBRARIES})
  SET(FEEL_LIBRARIES ${MadLib_LIBRARY} ${FEEL_LIBRARIES})
endif( MADLIB_FOUND )

#
# Eigen
#
set(EIGEN_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
add_subdirectory(contrib/eigen)
INCLUDE_DIRECTORIES( ${FEEL_SOURCE_DIR}/contrib/eigen )

#FIND_PACKAGE(Eigen2 REQUIRED)
#INCLUDE_DIRECTORIES( ${Eigen2_INCLUDE_DIR} )
#add_subdirectory(contrib/eigen)
#INCLUDE_DIRECTORIES( ${FEEL_SOURCE_DIR}/contrib/eigen )
#add_definitions( -DEIGEN_NO_STATIC_ASSERT )

#
# Scotch
#
#CheckIncludeFileCXX( ptscotch.h HAVE_PTSCOTCH_H )
#CheckIncludeFileCXX( scotch.h HAVE_SCOTCH_H )
#
# Metis
#
FIND_PACKAGE(Metis)
if ( METIS_FOUND )
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LINK_DIRECTORIES(${METIS_LIBRARIES})
  SET(FEEL_LIBRARIES ${METIS_LIBRARY} ${FEEL_LIBRARIES})
endif( METIS_FOUND )

#
# Petsc
#
FIND_PACKAGE( PETSc )
if ( PETSC_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
  SET(FEEL_LIBRARIES ${PETSC_LIBRARIES} ${FEEL_LIBRARIES})
  SET(BACKEND_PETSC petsc)
endif( PETSC_FOUND )

#
# parpack
#
FIND_LIBRARY(PARPACK_LIBRARY NAMES parpack)
if (PARPACK_LIBRARY)
    SET(PARPACK_LIBRARIES ${PARPACK_LIBRARY})
    SET(FEEL_LIBRARIES ${PARPACK_LIBRARIES} ${FEEL_LIBRARIES})
endif()
MARK_AS_ADVANCED( PARPACK_LIBRARY )


#
# SLEPc
#
FIND_PACKAGE( SLEPc )
if ( SLEPC_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${SLEPC_INCLUDES};${CMAKE_REQUIRED_INCLUDES}")
  INCLUDE_DIRECTORIES( ${SLEPC_INCLUDE_DIR} )
  SET(FEEL_LIBRARIES ${SLEPC_LIBRARIES} ${FEEL_LIBRARIES})
endif(SLEPC_FOUND)


#
# Trilinos
#
FIND_PACKAGE(Trilinos)
if ( TRILINOS_FOUND )
  INCLUDE_DIRECTORIES(${TRILINOS_INCLUDE_DIR})
  SET(FEEL_LIBRARIES ${TRILINOS_LIBRARIES} ${FEEL_LIBRARIES})
  SET(BACKEND_TRILINOS trilinos)
endif( TRILINOS_FOUND )

#
# OpenTURNS
#
FIND_PACKAGE( OpenTURNS )
if ( OPENTURNS_FOUND )
  MESSAGE(STATUS "OpenTURNS Libraries: ${OpenTURNS_LIBRARIES}")
  MESSAGE(STATUS "OpenTURNS Headers: ${OpenTURNS_INCLUDE_DIRS}")
  INCLUDE_DIRECTORIES(${OpenTURNS_INCLUDE_DIRS})
  SET(FEEL_LIBRARIES ${OpenTURNS_LIBRARIES} ${FEEL_LIBRARIES})
endif( OPENTURNS_FOUND )

#
# VTK
#
FIND_PACKAGE(VTK)
if ( VTK_FOUND )
  set(HAVE_VTK 1)
  SET(VTK_LIBRARIES "-lvtkRendering -lvtkGraphics -lvtkImaging -lvtkCommon" )
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
  MARK_AS_ADVANCED( VTK_DIR )
  SET(FEEL_LIBRARIES ${VTK_LIBRARIES} ${FEEL_LIBRARIES})
endif()

#
# Octave
#
FIND_PACKAGE(Octave)
if ( OCTAVE_FOUND )
  INCLUDE_DIRECTORIES( ${Octave_INCLUDE_DIR} )
endif( OCTAVE_FOUND )

#
#
#
INCLUDE_DIRECTORIES (
  ${FEEL_BUILD_DIR}/
  ${FEEL_SOURCE_DIR}/
  ${FEEL_SOURCE_DIR}/contrib/gmm/include

  ${PETSC_INCLUDE_DIR}
  ${PETSCCONF_INCLUDE_DIR}

  ${MPI_INCLUDE_PATH}
  ${BOOST_INCLUDE_PATH}
  ${LIBXML2_INCLUDE_DIR}
  )

LINK_DIRECTORIES(
  ${VTK_LIBRARY_DIRS}
  ${BOOST_LIBRARY_PATH}
  ${MPI_LIBRARY_PATH}
  )

#
# programs
#
find_program( GMSH gmsh DOC "Gmsh mesh generator" )
if ( GMSH )
  ADD_DEFINITIONS( -DHAVE_GMSH=1 -D_HAVE_GMSH_)
endif()
mark_as_advanced( GMSH )

# ${FEEL_SOURCE_DIR}/feel/feelcore
# ${FEEL_SOURCE_DIR}/feel/feelalg
# ${FEEL_SOURCE_DIR}/feel/feelmesh
# ${FEEL_SOURCE_DIR}/feel/feelpoly
# ${FEEL_SOURCE_DIR}/feel/feeldiscr
# ${FEEL_SOURCE_DIR}/feel/feelvf
# ${FEEL_SOURCE_DIR}/feel/feelfilters
# ${FEEL_SOURCE_DIR}/feel/feelsystem

#
# Python
#
FIND_PACKAGE(PythonInterp 2.2 REQUIRED)

#
# Feel
#
ADD_SUBDIRECTORY ( feel )

SET(FEEL_LIBRARIES
  feel++
  ${FEEL_LIBRARIES}
  )



SET(FEEL_MESH_MAX_ORDER "5" CACHE STRING "maximum geometrical order in templates to instantiate" )

OPTION(FEEL_ENABLE_DOCUMENTATION "enable feel documentation" OFF)
OPTION(FEEL_ENABLE_BENCHMARKS "enable feel benchmarks" OFF)
OPTION(FEEL_ENABLE_EXAMPLES "enable feel examples" OFF)
OPTION(FEEL_ENABLE_APPLICATIONS "enable feel applications" OFF)
OPTION(FEEL_ENABLE_TESTS "enable feel tests" OFF)
OPTION(FEEL_ENABLE_ALL "enable all feel modules" OFF)

if ( FEEL_ENABLE_ALL )
  set( FEEL_ENABLE_DOCUMENTATION ON )
  set( FEEL_ENABLE_BENCHMARKS ON )
  set( FEEL_ENABLE_EXAMPLES ON )
  set( FEEL_ENABLE_APPLICATIONS ON )
  set( FEEL_ENABLE_TESTS ON )
endif( FEEL_ENABLE_ALL )

OPTION(FEEL_MINIMAL_CONFIGURATION "enable feel minimal configuration" OFF)
IF( FEEL_MINIMAL_CONFIGURATION )
  set( FEEL_ENABLE_DOCUMENTATION OFF )
  set( BUILD_TESTING OFF )
  set( FEEL_ENABLE_DOCUMENTATION ON )
  set( ENABLE_INSTANTIATION_MODE OFF )
  UNSET( FEEL_INSTANTIATION_MODE CACHE )
  SET(FEEL_MESH_MAX_ORDER "1" CACHE STRING "maximum geometrical order in templates to instantiate" FORCE )
ENDIF( FEEL_MINIMAL_CONFIGURATION )

#
# Enable testing
#
INCLUDE(CTest)
ENABLE_TESTING()

if ( FEEL_ENABLE_DOCUMENTATION )
  ADD_SUBDIRECTORY ( doc )
endif ( FEEL_ENABLE_DOCUMENTATION )

if ( FEEL_ENABLE_BENCHMARKS )
  ADD_SUBDIRECTORY ( benchmarks )
endif()


if ( FEEL_ENABLE_EXAMPLES )
 ADD_SUBDIRECTORY ( examples )
endif()

if ( FEEL_ENABLE_APPLICATIONS )
 ADD_SUBDIRECTORY ( applications )
 if ( EXISTS ${FEEL_SOURCE_DIR}/research )
  ADD_SUBDIRECTORY ( research )
 endif()
endif( FEEL_ENABLE_APPLICATIONS )




IF(FEEL_ENABLE_TESTS )

  macro(feel_add_test testname)
    set(targetname test_${testname})
    set(filename test_${testname}.cpp)
    add_executable(${targetname} ${filename})
    target_link_libraries(${targetname} ${FEEL_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
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
  endmacro(feel_add_test)

  add_subdirectory( testsuite )

endif()

# generate configuration header
set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEEL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(FEEL_DATADIR ${CMAKE_INSTALL_PREFIX}/share/feel )
CONFIGURE_FILE(feelconfig.h.in feelconfig.h  @ONLY)
CONFIGURE_FILE(feelinfo.h.in feelinfo.h  @ONLY)

#
# Packaging
#
INCLUDE(InstallRequiredSystemLibraries)


foreach(includedir feelcore feelalg feelmesh feelpoly feelfilters feeldiscr feelvf feelmaterial feelsystem )
  FILE(GLOB files "feel/${includedir}/*.hpp")
  if ( NOT ENABLE_INSTANTIATION_MODE )
    FILE(GLOB cppfiles "feel/${includedir}/*.cpp")
  endif()
  INSTALL(FILES ${files} ${cppfiles} DESTINATION include/feel/${includedir} COMPONENT Devel)
endforeach()
FILE(GLOB files "feel/*.hpp")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
FILE(GLOB files "contrib/gmm/include/*.h")
INSTALL(FILES ${files} DESTINATION include/feel COMPONENT Devel)
FILE(GLOB files ${CMAKE_CURRENT_BINARY_DIR}/feelconfig.h)
INSTALL(FILES ${files} DESTINATION include/ COMPONENT Devel)
FILE(GLOB libso "${CMAKE_CURRENT_BINARY_DIR}/feel/libfeel++.so*")
INSTALL(FILES ${libso} DESTINATION lib/ COMPONENT Libraries)
foreach( programs_dir examples/generic examples/heat examples/solid doc/tutorial applications/polyvis benchmarks/convergence )
  FILE(GLOB programs "${CMAKE_CURRENT_BINARY_DIR}/${programs_dir}/feel_*")
  INSTALL(FILES ${programs} DESTINATION bin/ COMPONENT Applications)
endforeach()
# documentation and examples
FILE(GLOB pdfs "${CMAKE_CURRENT_BINARY_DIR}/doc/tutorial/feel-manual.pdf")
INSTALL(FILES ${pdfs} DESTINATION share/doc/feel/ COMPONENT Documentation)
FILE(GLOB examples "${CMAKE_CURRENT_SOURCE_DIR}/doc/tutorial/*.*pp")
INSTALL(FILES ${examples} DESTINATION share/doc/feel/examples/ COMPONENT Examples)


##
## Archive generation using cpack
##
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Feel++")
SET(CPACK_PACKAGE_VENDOR "Christophe Prud'homme")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FEEL_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FEEL_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FEEL_VERSION_MICRO}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "feel")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "CPackSourceConfig.cmake")
if ( EXISTS ${FEEL_SOURCE_DIR}/applications/opus )
  #SET(CPACK_SOURCE_PACKAGE_FILE_NAME "feel+opus-${FEEL_VERSION_MAJOR}.${FEEL_VERSION_MINOR}.${FEEL_VERSION_MICRO}${FEEL_VERSION_EXTRA}")
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "feel++-${FEEL_VERSION_MAJOR}.${FEEL_VERSION_MINOR}.${FEEL_VERSION_MICRO}${FEEL_VERSION_EXTRA}")
else()
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "feel++-${FEEL_VERSION_MAJOR}.${FEEL_VERSION_MINOR}.${FEEL_VERSION_MICRO}${FEEL_VERSION_EXTRA}")
endif()

SET(CPACK_SOURCE_STRIP_FILES "")
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git/;\\\\.gitignore;/\\\\.svn;"
  "/admin/;/tools/;/Templates/;"
  "/auto/;/auto.*/;"
  "/TAGS;/#.*;/.*~;/.cvsignore;/.bzrignore"
  "${PROJECT_SOURCE_DIR}/benchmarks/turek/"
  "${PROJECT_SOURCE_DIR}/benchmarks/stokes/"
  "${PROJECT_SOURCE_DIR}/benchmarks/ethiersteinman/"
  "${PROJECT_SOURCE_DIR}/benchmarks/kovasznay/"
  "${PROJECT_SOURCE_DIR}/doc/poster/"
  "${PROJECT_SOURCE_DIR}/doc/tutorial/pdfs/"
  "${PROJECT_SOURCE_DIR}/doc/figures/backgrounds/"
  "${PROJECT_SOURCE_DIR}/doc/figures/logos/"
  "${PROJECT_SOURCE_DIR}/examples/fluid/"
  "${PROJECT_SOURCE_DIR}/examples/levelset/"
  "${PROJECT_SOURCE_DIR}/examples/pbeq/"
  "${PROJECT_SOURCE_DIR}/research/"  )

#if ( NOT EXISTS ${FEEL_SOURCE_DIR}/applications/opus/ )
if ( EXISTS ${FEEL_SOURCE_DIR}/applications/opus/ )
  set(CPACK_SOURCE_IGNORE_FILES "${PROJECT_SOURCE_DIR}/applications/opus/;${CPACK_SOURCE_IGNORE_FILES}" )
endif()

include( CPack )



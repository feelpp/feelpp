###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-05-27
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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
# GFLAGS_FOUND - system has GFLAGS
# GFLAGS_INCLUDE_DIR - headers location
# GFLAGS_LIBRARIES - libraries

# try to find gflags headers
# and set GFLAGS_INCLUDE_DIR and GFLAGS_LIBRARIES

# if the gflags source directory exists then we use it even if an elligible
# version of gflags is available on the system
if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gflags )
  # try local version if system version is not present
  FIND_PATH(GFLAGS_INCLUDE_DIR gflags/gflags.h
    ${CMAKE_BINARY_DIR}/contrib/gflags/include
    NO_DEFAULT_PATH
    )
  message(STATUS "Gflags first pass: ${GFLAGS_INCLUDE_DIR}")

  if (NOT GFLAGS_INCLUDE_DIR )

    execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/contrib/gflags-compile)
    if(${CMAKE_SOURCE_DIR}/contrib/gflags/configure.ac IS_NEWER_THAN ${CMAKE_BINARY_DIR}/contrib/gflags-compile/configure)
      message(STATUS "Building gflags in ${CMAKE_BINARY_DIR}/contrib/gflags-compile...")
      if (FEELPP_USE_STATIC_LINKAGE )
        message(STATUS "GFlags: use static linkage")
        execute_process(
          COMMAND ${FEELPP_HOME_DIR}/contrib/gflags/configure --prefix=${CMAKE_BINARY_DIR}/contrib/gflags  --enable-static --disable-shared  CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gflags-compile
          #      OUTPUT_QUIET
          OUTPUT_FILE "gflags-configure"
          )
      else()
        execute_process(
          COMMAND ${FEELPP_HOME_DIR}/contrib/gflags/configure --prefix=${CMAKE_BINARY_DIR}/contrib/gflags  LDFLAGS=-dynamic CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gflags-compile
          #      OUTPUT_QUIET
          OUTPUT_FILE "gflags-configure"
          )
      endif(FEELPP_USE_STATIC_LINKAGE)
    endif()

    set(GFLAGS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/contrib/gflags/include)

  endif()

  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gflags/ )
    if(${CMAKE_SOURCE_DIR}/contrib/gflags/src/gflags/gflags.h IS_NEWER_THAN ${CMAKE_BINARY_DIR}/contrib/gflags/include/gflags/gflags.h)
      message(STATUS "Installing gflags in ${CMAKE_BINARY_DIR}/contrib/gflags...")
      if ( APPLE AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"  OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang"))
        execute_process(
          COMMAND make -k -j${NProcs2} install CXXFLAGS=-stdlib=libc++
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gflags-compile
          #  OUTPUT_QUIET
          OUTPUT_FILE "gflags-install"
          )
      else()
        execute_process(
          COMMAND make -k -j${NProcs2} install
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/gflags-compile
          #  OUTPUT_QUIET
          OUTPUT_FILE "gflags-install"
          )
      endif()
    endif()
  endif()


  FIND_LIBRARY(GFLAGS_LIBRARY
    NAMES feelpp_gflags
    PATHS
    ${CMAKE_BINARY_DIR}/contrib/gflags/lib64/
    ${CMAKE_BINARY_DIR}/contrib/gflags/lib/
    NO_DEFAULT_PATH
    )
else( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gflags )

  FIND_PATH(GFLAGS_INCLUDE_DIR gflags/gflags.h
    $ENV{FEELPP_DIR}/include/feel
    NO_DEFAULT_PATH)

  FIND_PATH(GFLAGS_INCLUDE_DIR gflags/gflags.h
    $ENV{FEELPP_DIR}/include/feel
    /usr/include/feel
    /usr/local/include/feel
    /opt/local/include/feel
    NO_DEFAULT_PATH)
  message(STATUS "Gflags/system: ${GFLAGS_INCLUDE_DIR}")
  FIND_LIBRARY(GFLAGS_LIBRARY  NAMES feelpp_gflags  PATHS   $ENV{FEELPP_DIR}/lib  NO_DEFAULT_PATH)
  FIND_LIBRARY(GFLAGS_LIBRARY  NAMES feelpp_gflags    )
endif( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gflags )

string(REPLACE "include" "" GFLAGS_DIR ${GFLAGS_INCLUDE_DIR} )

set(GFLAGS_LIBRARIES ${GFLAGS_LIBRARY})
message(STATUS "Gflags includes: ${GFLAGS_INCLUDE_DIR} Libraries: ${GFLAGS_LIBRARIES} Dir: ${GFLAGS_DIR}" )


# handle the QUIETLY and REQUIRED arguments and set GFLAGS_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GFLAGS DEFAULT_MSG GFLAGS_INCLUDE_DIR GFLAGS_LIBRARIES GFLAGS_DIR )

mark_as_advanced (GFLAGS_INCLUDE_DIR GFLAGS_LIBRARIES GFLAGS_DIR GFLAGS_LIBRARY)

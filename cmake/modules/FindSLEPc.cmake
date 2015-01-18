# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-01-22
#
#  Copyright (C) 2010 Université Joseph Fourier
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

set(DARWIN_FLAVORS real complex)

IF ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
  set( DARWIN_FLAVORS darwin-cxx-debug arch-darwin-cxx-debug arch-darwin-cxx-opt darwin-cxx-opt   ${DARWIN_FLAVORS})
ELSE()
  set( DARWIN_FLAVORS darwin-cxx-opt  arch-darwin-cxx-opt darwin-cxx-debug arch-darwin-cxx-debug  ${DARWIN_FLAVORS})
ENDIF()

foreach( debian_arches linux kfreebsd )
  IF ( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
    set( DEBIAN_FLAVORS ${debian_arches}-gnu-c-debug ${debian_arches}-gnu-c-opt ${DEBIAN_FLAVORS})
  ELSE()
    set( DEBIAN_FLAVORS ${debian_arches}-gnu-c-opt ${debian_arches}-gnu-c-debug ${DEBIAN_FLAVORS})
  ENDIF()
endforeach()

set(PETSC_VERSIONS 3.5.2 3.5.1 3.5.0 3.4.4 3.4.3 3.4.2 3.3 3.2 )

if ( NOT SLEPC_DIR )
  foreach( version ${PETSC_VERSIONS} )
    foreach ( flavor ${DEBIAN_FLAVORS} ${DARWIN_FLAVORS})
      #message(STATUS "checking version ${version} for file ${flavor}/include/petsc.h...")
      find_path (SLEPC_DIR include/slepc.h
        PATHS
        /usr/lib/slepcdir/${version}/${flavor}
        /usr/local/Cellar/slepc/${version}/${flavor}
        NO_DEFAULT_PATH
        DOC "SLEPc Directory")
    endforeach()
  endforeach()
endif()
message(STATUS "SLEPc Dir: ${SLEPC_DIR}")


find_path (SLEPC_DIR include/slepc.h
  HINTS ENV SLEPC_DIR
  PATHS
  /usr/lib/slepc
  /usr/lib/slepcdir/3.4.4 # Debian
  /usr/lib/slepcdir/3.4.3 # Debian
  /usr/lib/slepcdir/3.4.2 # Debian
  /usr/lib/slepcdir/3.2 # Debian
  /usr/lib/slepcdir/3.1 # Debian
  /usr/lib/slepcdir/3.0.0 # Debian
  /opt/local/lib/petsc # macports
  /opt/local/lib/slepc # macports
  # Homebrew
  /usr/local/Cellar/slepc/3.5.2
  /usr/local/Cellar/slepc/3.5.1
  /usr/local/Cellar/slepc/3.5.0
  /usr/local/Cellar/slepc/3.4.4
  /usr/local/Cellar/slepc/3.4.3
  $ENV{HOME}/slepc
  DOC "SLEPc Directory")



SET(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include/")
CHECK_INCLUDE_FILE( ${SLEPC_INCLUDE_DIR}/slepc.h FEELPP_HAS_SLEPC_H )

if (SLEPC_DIR AND NOT PETSC_ARCH)
  set (_slepc_arches
    $ENV{PETSC_ARCH}                   # If set, use environment variable first
    ${DEBIAN_FLAVORS}  # Debian defaults
    ${DARWIN_FLAVORS}  # Darwin defaults
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
  set (slepcconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_slepc_arches})
    if (NOT PETSC_ARCH)
      find_path (slepcconf slepcconf.h
	HINTS ${SLEPC_DIR}
	PATH_SUFFIXES ${arch}/include bmake/${arch}
	NO_DEFAULT_PATH)
      if (slepcconf)
	set (PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture")
      endif (slepcconf)
    endif (NOT PETSC_ARCH)
  endforeach (arch)
  set (slepcconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif (SLEPC_DIR AND NOT PETSC_ARCH)
message(STATUS "SLEPc arch: ${PETSC_ARCH}")

#CHECK_INCLUDE_FILE( ${SLEPC_DIR}/$ENV{PETSC_ARCH}/include/slepcconf.h FEELPP_HAS_SLEPCCONF_H )
#if (FEELPP_HAS_SLEPCCONF_H)
  SET(SLEPC_INCLUDE_DIR ${SLEPC_DIR}/${PETSC_ARCH}/include ${SLEPC_INCLUDE_DIR})
#endif()
FIND_LIBRARY(SLEPC_LIB_SLEPC slepc
  HINTS
  ${SLEPC_DIR}/${PETSC_ARCH}/lib
  ${SLEPC_DIR}/lib )
SET(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC} CACHE STRING "SLEPc libraries" FORCE)
# handle the QUIETLY and REQUIRED arguments and set Madlib_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (SLEPC DEFAULT_MSG
  SLEPC_DIR SLEPC_LIBRARIES
  )

if ( SLEPC_FOUND )
  MESSAGE( STATUS "SLepc found: ${SLEPC_DIR}" )
  set(FEELPP_HAS_SLEPC 1)
  set(SLEPC_INCLUDES ${SLEPC_INCLUDE_DIR} CACHE STRING "SLEPc include path" FORCE)
endif()

MARK_AS_ADVANCED( SLEPC_DIR SLEPC_LIB_SLEPC SLEPC_INCLUDES SLEPC_LIBRARIES )

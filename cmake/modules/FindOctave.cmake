# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-04-10
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

include (CheckIncludeFileCXX)
include (FindPackageHandleStandardArgs)
find_package(Readline REQUIRED)

# test if variables are not already in cache
if (NOT (ENABLE_Octave
      AND Octave_INCLUDE_DIR
      AND Octave_INCLUDE_DIRS
      AND Octave_LIBRARY
      AND Octave_LIBRARIES))

  find_program(OCTAVE octave)
  if (OCTAVE)
    message(STATUS "Found Octave: ${OCTAVE}")
    set(ENABLE_Octave ON CACHE BOOL "Enable Octave bindings" FORCE)
  else(OCTAVE)
    message(STATUS "WARNING: "
      "octave not found. Disabling octave bindings")
    set(ENABLE_Octave OFF CACHE BOOL "Enable Octave bindings" FORCE)
  endif(OCTAVE)

  if(ENABLE_Octave)
    #OCTAVE_VERSION is the (dotted triplet) octave version.
    execute_process(
      COMMAND ${OCTAVE} --version
      OUTPUT_VARIABLE _OCTAVE_VERSION
      )
    string(REGEX REPLACE
      "^.*version ([0-9]\\.[0-9]\\.[0-9]*).*$"
      "\\1"
      OCTAVE_VERSION
      ${_OCTAVE_VERSION}
      )
    message(STATUS "OCTAVE_VERSION: ${OCTAVE_VERSION}")
    # Logic that depends on octave version
    # transform_version(NUMERICAL_OCTAVE_TESTING_MINIMUM_VERSION "3.2.0")
    # transform_version(NUMERICAL_OCTAVE_VERSION "${OCTAVE_VERSION}")
    set(NUMERICAL_OCTAVE_TESTING_MINIMUM_VERSION "3.2.0")
    set(NUMERICAL_OCTAVE_VERSION "${OCTAVE_VERSION}")
    if(
        NUMERICAL_OCTAVE_VERSION
        LESS
        "${NUMERICAL_OCTAVE_TESTING_MINIMUM_VERSION}"
        )
      message(STATUS "WARNING: "
        "Opus require octave version 3.2 or greater. Disabling octave bindings")
      set(ENABLE_Octave OFF CACHE BOOL "Enable Octave bindings" FORCE)
    endif(
      NUMERICAL_OCTAVE_VERSION
      LESS
      "${NUMERICAL_OCTAVE_TESTING_MINIMUM_VERSION}"
      )
  endif(ENABLE_Octave)
  MESSAGE(STATUS "Octave version: " ${NUMERICAL_OCTAVE_VERSION} ${OCTAVE_VERSION})
  # set include dir
  if (NOT Octave_INCLUDE_DIR)
    find_path (Octave_INCLUDE_DIR
      NAMES
      oct.h
      PATHS
      /usr/include
      /usr/local/include
      /opt/local/include
      /sw/include
      PATH_SUFFIXES
      octave
      octave-${OCTAVE_VERSION}
      octave-${OCTAVE_VERSION}/octave
      DOC
      "Octave include directory"
      )
  endif ()

  # dependencies includes
  if (NOT Octave_INCLUDE_DIRS)
    set (Octave_INCLUDE_DIRS ${Octave_INCLUDE_DIR})
    #  list (APPEND Octave_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
  endif ()

  if (NOT Octave_LIBRARIES)
    find_library(
      Octave_LIBRARY
      octave
      PATH_SUFFIXES octave-${OCTAVE_VERSION}
      )

    find_library(
      Octinterp_LIBRARY
      octinterp
      PATH_SUFFIXES octave-${OCTAVE_VERSION}
      )

    find_library(
      CRUFT_LIBRARY
      cruft
      PATH_SUFFIXES octave-${OCTAVE_VERSION}
      )
    MESSAGE(STATUS "Octave library: ${OCTAVE_LIBRARY}")
    MESSAGE(STATUS "Octinterp library: ${Octinterp_LIBRARY}")
    MESSAGE(STATUS "Cruft library: ${CRUFT_LIBRARY}")

    # find dependent libraries
    set (Octave_LIBRARIES ${Octave_LIBRARY} ${Octinterp_LIBRARY})
    list (APPEND Octave_LIBRARIES ${READLINE_LIBRARY})
    list (APPEND Octave_LIBRARIES ${CRUFT_LIBRARY})
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set Octave_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args (Octave DEFAULT_MSG
  Octave_LIBRARY
  Octave_INCLUDE_DIR
  Octave_INCLUDE_DIRS
  Octave_LIBRARIES
  )

mark_as_advanced (
  Octave_LIBRARY
  Octave_INCLUDE_DIR
  Octave_INCLUDE_DIRS
  Octave_LIBRARIES
  )

if ( OCTAVE_FOUND )
  SET( FEELPP_HAS_OCT_H 1 )
endif()

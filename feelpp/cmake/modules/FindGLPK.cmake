# -*- mode: cmake -*-
#
#  This file is part of the Feel++ library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2010-02-10
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
# this files defines
#  - GLPK_INCLUDE_DIR
#  - GLPK_LIBRARIES
#  - GLPK_FOUND

INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(glpk.h FEELPP_HAS_GLPK_H)

if(DEFINED ENV{SPACK_ENV} AND NOT "$ENV{SPACK_ENV}" STREQUAL "")
  list(APPEND _GLPK_PREFIX_HINTS ${CMAKE_PREFIX_PATH})
  if(DEFINED ENV{CMAKE_PREFIX_PATH} AND NOT "$ENV{CMAKE_PREFIX_PATH}" STREQUAL "")
    set(_GLPK_ENV_PREFIX_HINTS "$ENV{CMAKE_PREFIX_PATH}")
    if(NOT WIN32)
      string(REPLACE ":" ";" _GLPK_ENV_PREFIX_HINTS "${_GLPK_ENV_PREFIX_HINTS}")
    endif()
    list(APPEND _GLPK_PREFIX_HINTS ${_GLPK_ENV_PREFIX_HINTS})
  endif()
  list(FILTER _GLPK_PREFIX_HINTS EXCLUDE REGEX "^$")
  list(REMOVE_DUPLICATES _GLPK_PREFIX_HINTS)

  set(_GLPK_REAL_PREFIX_HINTS)
  foreach(_GLPK_PREFIX_HINT IN LISTS _GLPK_PREFIX_HINTS)
    if(EXISTS "${_GLPK_PREFIX_HINT}")
      get_filename_component(_GLPK_REAL_PREFIX_HINT "${_GLPK_PREFIX_HINT}" REALPATH)
      list(APPEND _GLPK_REAL_PREFIX_HINTS "${_GLPK_REAL_PREFIX_HINT}")
    endif()
  endforeach()

  foreach(_GLPK_CACHE_VAR GLPK_LIB GLPK_INCLUDE_DIR)
    if(DEFINED ${_GLPK_CACHE_VAR} AND NOT "${${_GLPK_CACHE_VAR}}" MATCHES "-NOTFOUND$")
      get_filename_component(_GLPK_REAL_CACHE_VALUE "${${_GLPK_CACHE_VAR}}" REALPATH)
      set(_GLPK_CACHE_IN_ACTIVE_PREFIX FALSE)
      foreach(_GLPK_REAL_PREFIX_HINT IN LISTS _GLPK_REAL_PREFIX_HINTS)
        if(_GLPK_REAL_CACHE_VALUE MATCHES "^${_GLPK_REAL_PREFIX_HINT}(/|$)")
          set(_GLPK_CACHE_IN_ACTIVE_PREFIX TRUE)
        endif()
      endforeach()
      if(NOT _GLPK_CACHE_IN_ACTIVE_PREFIX)
        message(STATUS "[glpk] ignoring cached ${_GLPK_CACHE_VAR}='${${_GLPK_CACHE_VAR}}' outside active Spack prefixes")
        unset(${_GLPK_CACHE_VAR} CACHE)
        unset(${_GLPK_CACHE_VAR})
      endif()
    endif()
  endforeach()
endif()

FIND_LIBRARY( GLPK_LIB glpk
  PATHS ${_GLPK_PREFIX_HINTS} $ENV{GLPK_DIR} /usr
  PATH_SUFFIXES lib lib64)
SET(GLPK_LIBRARIES ${GLPK_LIB} )

FIND_PATH(GLPK_INCLUDE_DIR
  glpk.h
  PATHS ${_GLPK_PREFIX_HINTS} $ENV{GLPK_DIR} /usr
  PATH_SUFFIXES include include/glpk
  DOC "Directory where GLPK header files are stored" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK "Could not find GLPK " GLPK_INCLUDE_DIR GLPK_LIBRARIES)
# show the BERKELEY_DB_INCLUDE_DIR and BERKELEY_DB_LIBRARIES variables only in the advanced view
MARK_AS_ADVANCED(GLPK_INCLUDE_DIR GLPK_LIBRARIES )

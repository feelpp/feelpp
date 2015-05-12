###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
#       Date: 2014-12-18
#
#  Copyright (C) 2014-2015 Feel++ Consortium
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

#
# METIS
#
OPTION( FEELPP_ENABLE_METIS "Enable METIS" OFF )

if ( FEELPP_ENABLE_METIS )
  #
  # Metis
  #
  FIND_PACKAGE(Metis)
  if ( METIS_FOUND )
    INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
    #  LINK_DIRECTORIES(${METIS_LIBRARIES})
    SET(FEELPP_LIBRARIES ${METIS_LIBRARY} ${FEELPP_LIBRARIES})
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Metis" )
  endif( METIS_FOUND )
  
  # metis
  FIND_LIBRARY(METIS_LIBRARY
    NAMES
    metis
    PATHS
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    #    "/opt/local/lib"
    )
  message(STATUS "[feelpp] Metis: ${METIS_LIBRARY}" )
  IF( METIS_LIBRARY )
    SET(FEELPP_LIBRARIES ${METIS_LIBRARY} ${FEELPP_LIBRARIES})
  ENDIF()
  
  
  if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/metis/ )
    message( FATAL_ERROR "Please make sure that git submodule contrib/metis is available. Run `git submodule update --init --recursive contrib/metis`")
  endif()

  FIND_PATH(METIS_INCLUDE_DIR metis.h HINTS ${FEELPP_SOURCE_DIR}/contrib $ENV{METIS_DIR} ${METIS_INCLUDE_DIR} PATH_SUFFIXES metis/include)
  if( METIS_INCLUDE_DIR )
    INCLUDE_DIRECTORIES( ${METIS_INCLUDE_DIR} )
    SET(FEELPP_HAS_METIS 1)
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} METIS" )
    ADD_DEFINITIONS( -DFEELPP_HAS_METIS )
    add_subdirectory(contrib/metis)
  endif()

endif()

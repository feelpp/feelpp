###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-01-28
#
#  Copyright (C) 2014 Feel++ Consortium
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

find_path (SCOTCH_DIR include/ptscotch.h
  HINTS ENV PTSCOTCH_DIR
  PATHS
  $ENV{PETSC_DIR}/
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/
  /opt/local/lib/petsc/
  /usr/local/opt/scotch5/
  )

SET(SCOTCH_INCLUDE_DIR "${SCOTCH_DIR}/include/")
CHECK_INCLUDE_FILE( ${SCOTCH_INCLUDE_DIR}/scotch.h FEELPP_HAS_SCOTCH_H )

FIND_LIBRARY(PTSCOTCHERREXIT_LIBRARY
  NAMES
  ptscotcherrexit
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  /opt/local/lib/petsc/lib
  $ENV{PTSCOTCH_DIR}/lib
  /usr/local/opt/scotch5/lib
  )
message(STATUS "PTScotcherrexit: ${PTSCOTCHERREXIT_LIBRARY}" )
IF( PTSCOTCHERREXIT_LIBRARY )
  SET(SCOTCH_LIBRARIES ${PTSCOTCHERREXIT_LIBRARY} ${SCOTCH_LIBRARIES})
ENDIF()

FIND_LIBRARY(PTSCOTCHERR_LIBRARY
  NAMES
  ptscotcherr
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  /opt/local/lib/petsc/lib
  $ENV{PTSCOTCH_DIR}/lib
  /usr/local/opt/scotch5/lib
  )
message(STATUS "PTScotcherr: ${PTSCOTCHERR_LIBRARY}" )
IF( PTSCOTCHERR_LIBRARY )
  #message(STATUS "PTScotch: ${PTSCOTCH_LIBRARY}" )
  SET(SCOTCH_LIBRARIES ${PTSCOTCHERR_LIBRARY} ${SCOTCH_LIBRARIES})
ENDIF()

FIND_LIBRARY(PTSCOTCH_LIBRARY
  NAMES
  ptscotch
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  /opt/local/lib/petsc/lib
  $ENV{PTSCOTCH_DIR}/lib
  /usr/local/opt/scotch5/lib
  )
message(STATUS "PTScotch: ${PTSCOTCH_LIBRARY}" )

IF( PTSCOTCH_LIBRARY )
  SET(SCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${SCOTCH_LIBRARIES})
ENDIF()

FIND_LIBRARY(PTESMUMPS_LIBRARY
  NAMES
  ptesmumps
  PATHS
  $ENV{PETSC_DIR}/lib
  $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
  /opt/local/lib/petsc/lib
  $ENV{PTSCOTCH_DIR}/lib
  /usr/local/opt/scotch5/lib
  )
message(STATUS "PTESMUMPS: ${PTESMUMPS_LIBRARY}" )
IF( PTESMUMPS_LIBRARY )
  set( FEELPP_HAS_MUMPS 1 )
  SET(SCOTCH_LIBRARIES ${PTESMUMPS_LIBRARY} ${SCOTCH_LIBRARIES})
ENDIF()

# handle the QUIETLY and REQUIRED arguments and set Madlib_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (SCOTCH DEFAULT_MSG
  SCOTCH_DIR SCOTCH_LIBRARIES
  )

if ( SCOTCH_FOUND )
  MESSAGE( STATUS "Scotch found: ${SCOTCH_DIR}" )
  set(FEELPP_HAS_SCOTCH 1)
  set(SCOTCH_INCLUDES ${SCOTCH_INCLUDE_DIR} CACHE STRING "Scotch include path" FORCE)
endif()

MARK_AS_ADVANCED( SCOTCH_DIR SCOTCH_INCLUDES SCOTCH_LIBRARIES )

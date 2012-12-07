###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-09-07
#
#  Copyright (C) 2012 Université de Strasbourg
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
FIND_LIBRARY(MUMPS_COMMON_LIBRARY
    NAMES
    mumps_common
    PATHS
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    $ENV{MUMPS_DIR}/lib
)


FIND_LIBRARY(DMUMPS_LIBRARY
    NAMES
    dmumps
    PATHS
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    $ENV{MUMPS_DIR}/lib
)

message(STATUS "Mumps: ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY}" )
if ( MUMPS_COMMON_LIBRARY AND DMUMPS_LIBRARY )
  SET(MUMPS_LIBRARIES ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY})
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS "Could not find MUMPS " MUMPS_LIBRARIES)
MARK_AS_ADVANCED(MUMPS_LIBRARIES)

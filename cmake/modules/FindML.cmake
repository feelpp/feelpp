###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-16
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
include (FindPackageHandleStandardArgs)

option(FEELPP_ENABLE_ML_LIBRARY "Enables ML library in Feel++" ON )
if ( FEELPP_ENABLE_ML_LIBRARY )
  INCLUDE(CheckIncludeFileCXX)

  FIND_LIBRARY(ML_LIBRARY
    NAMES
    ml
    PATHS
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    ${PETSC_DIR}/lib
    NO_DEFAULT_PATH
    )

  FIND_LIBRARY(ML_LIBRARY
    NAMES
    ml
    PATHS
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    ${PETSC_DIR}/lib
    /opt/local/lib/petsc/lib
    )
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS (ML DEFAULT_MSG  ML_INCLUDE_DIR ML_LIBRARIES )

mark_as_advanced( ML_INCLUDE_DIR )
mark_as_advanced( ML_LIBRARY )

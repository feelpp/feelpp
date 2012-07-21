# -*- mode: cmake -*-
#
#  This file is part of the Feel++ library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
#  - MATHEVAL_INCLUDE_DIR
#  - MATHEVAL_LIBRARIES
#  - MATHEVAL_FOUND

INCLUDE(CheckIncludeFileCXX)


FIND_LIBRARY( MATHEVAL_LIB matheval PATHS /usr/lib /opt/local/lib  $ENV{MATHEVAL_DIR}/lib)
SET(MATHEVAL_LIBRARIES ${MATHEVAL_LIB} )

FIND_PATH(MATHEVAL_INCLUDE_DIR
  matheval.h
  PATHS /usr/include/ /usr/include/matheval /opt/local/include/matheval /usr/local/include/matheval  $ENV{MATHEVAL_DIR}/include/ANN
  DOC "Directory where MATHEVAL header files are stored" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MATHEVAL "Could not find MATHEVAL " MATHEVAL_INCLUDE_DIR MATHEVAL_LIBRARIES)
MARK_AS_ADVANCED(MATHEVAL_INCLUDE_DIR MATHEVAL_LIBRARIES )

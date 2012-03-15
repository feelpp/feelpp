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
#  - ANN_INCLUDE_DIR
#  - ANN_LIBRARIES
#  - ANN_FOUND

INCLUDE(CheckIncludeFileCXX)
# CHECK_INCLUDE_FILE_CXX(ANN/ANN.h FEELPP_HAS_ANN_H)


FIND_LIBRARY( ANN_LIB ann ANN PATHS /usr/lib /opt/local/lib  $ENV{ANN_DIR}/lib)
SET(ANN_LIBRARIES ${ANN_LIB} )

FIND_PATH(ANN_INCLUDE_DIR
  ANN.h
  PATHS /usr/include/ /usr/include/ANN /opt/local/include/ANN /usr/local/include/ANN  $ENV{ANN_DIR}/include/ANN
  DOC "Directory where ANN header files are stored" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ANN "Could not find ANN " ANN_INCLUDE_DIR ANN_LIBRARIES)
# show the BERKELEY_DB_INCLUDE_DIR and BERKELEY_DB_LIBRARIES variables only in the advanced view
MARK_AS_ADVANCED(ANN_INCLUDE_DIR ANN_LIBRARIES )

# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-03-15
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
include (FindPackageHandleStandardArgs)
INCLUDE(CheckIncludeFileCXX)

FIND_PATH(Eigen2_INCLUDE_DIR
  Eigen/Eigen
  PATHS /usr/include/ /usr/include/eigen2
  DOC "Directory where Eigen2 header files are stored" )


CHECK_INCLUDE_FILE_CXX(Eigen/Eigen FEELPP_HAS_EIGEN_EIGEN)
# handle the QUIETLY and REQUIRED arguments and set Eigen2_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (Eigen2 DEFAULT_MSG
  Eigen2_INCLUDE_DIR
  )

if ( EIGEN2_FOUND )
  MESSAGE( STATUS "Eigen2 found: ${Eigen2_INCLUDE_DIR}" )
  set( FEELPP_HAS_EIGEN_EIGEN 1 )
endif()

mark_as_advanced( Eigen2_INCLUDE_DIR )

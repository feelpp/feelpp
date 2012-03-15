# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-08-07
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

IF ( NOT TBB_INCLUDE_DIR)
  FIND_PATH(TBB_INCLUDE_DIR
    parallel_for.h
    PATH_SUFFIXES
    tbb
    PATHS /usr/include/ $ENV{TBB_INCLUDE_DIR}
    DOC "Directory where tbb header files are stored" )
ENDIF()


CHECK_INCLUDE_FILE_CXX(tbb.h FEELPP_HAS_TBB_H)

IF (NOT TBB_LIBRARIES)
  find_library(TBB_LIBRARY tbb PATHS $ENV{TBB_LIB_DIR})
  set (TBB_LIBRARIES ${TBB_LIBRARY})
ENDIF()

# handle the QUIETLY and REQUIRED arguments and set Eigen2_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (TBB DEFAULT_MSG
  TBB_INCLUDE_DIR
  TBB_LIBRARIES
  )

if ( TBB_FOUND )
  MESSAGE( STATUS "TBB found: headers: ${TBB_INCLUDE_DIR}, libs: ${TBB_LIBRARIES}" )
  SET( FEELPP_HAS_TBB 1 )
  SET( FEELPP_HAS_TBB_H 1 )
endif()

mark_as_advanced( TBB_INCLUDE_DIR )

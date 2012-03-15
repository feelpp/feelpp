# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-07-28
#
#  Copyright (C) 2010 Université de Grenoble 1 (Joseph Fourier)
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

FIND_PATH(MadLib_INCLUDE_DIR
  MAdLib.h
  PATHS /usr/include/ /usr/include/MAdLib/ /usr/include/madlib
  DOC "Directory where Madlib header files are stored" )


CHECK_INCLUDE_FILE_CXX(MadLib.h FEELPP_HAS_MADLIB_H)

FIND_LIBRARY(MadLib_LIBRARY MAdLib
  /usr/lib
  /usr/local/lib
  )
# handle the QUIETLY and REQUIRED arguments and set Madlib_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (MadLib DEFAULT_MSG
  MadLib_INCLUDE_DIR MadLib_LIBRARY
  )

if ( MadLib_FOUND )
  MESSAGE( STATUS "Madlib found: ${Madlib_INCLUDE_DIR}" )
endif()

mark_as_advanced( MadLib_INCLUDE_DIR )
mark_as_advanced( MadLib_LIBRARY )

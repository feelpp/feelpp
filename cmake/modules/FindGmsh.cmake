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

FIND_PATH(GMSH_INCLUDE_DIR
  Gmsh.h
  PATHS /usr/include/ /usr/include/gmsh/ /usr/local/include/gmsh /opt/local/include/gmsh $ENV{GMSH_DIR}/include/gmsh
  DOC "Directory where GMSH header files are stored" )


#CHECK_INCLUDE_FILE_CXX(Gmsh.h HAVE_GMSH_H)

FIND_LIBRARY(GMSH_LIBRARY NAMES Gmsh gmsh-2.5.1 gmsh
  PATH
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  $ENV{GMSH_DIR}/lib
  )

find_program( GMSH_EXECUTABLE gmsh DOC "GMSH mesh generator" )

# handle the QUIETLY and REQUIRED arguments and set Madlib_FOUND to TRUE if
# all listed variables are TRUE
FIND_PACKAGE_HANDLE_STANDARD_ARGS (GMSH DEFAULT_MSG
  GMSH_INCLUDE_DIR GMSH_LIBRARY GMSH_EXECUTABLE
  )

if ( GMSH_FOUND )
  MESSAGE( STATUS "GMSH found: header(${GMSH_INCLUDE_DIR}) lib(${GMSH_LIBRARY}) executable(${GMSH_EXECUTABLE})" )
endif()

mark_as_advanced( GMSH_INCLUDE_DIR )
mark_as_advanced( GMSH_LIBRARY )
mark_as_advanced( GMSH_EXECUTABLE )



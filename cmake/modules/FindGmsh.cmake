# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-07-28
#
#  Copyright (C) 2010 Universit� de Grenoble 1 (Joseph Fourier)
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

# some of the find_* commands are duplicated with the first instance having NO_DEFAULT_PATH
# this is to ensure that if GMSH_DIR is defined then cmake will pick the Gmsh version in GMSH_DIR
# otherwise in the second instance it will pich the standard version installed on the system if
# it is available

find_program( GMSH_EXECUTABLE
  NAMES gmsh
  PATHS
  $ENV{GMSH_DIR}/bin
  ${CMAKE_BINARY_DIR}/contrib/gmsh/bin
  PATH_SUFFIXES bin
  DOC "GMSH mesh generator" )

option(FEELPP_ENABLE_GMSH_LIBRARY "Enables Gmsh library in Feel++" ON )
if ( FEELPP_ENABLE_GMSH_LIBRARY )
  INCLUDE(CheckIncludeFileCXX)

  FIND_PATH(GMSH_INCLUDE_PATH
    Gmsh.h Context.h GModel.h
    HINTS
    $ENV{GMSH_DIR}
    ${CMAKE_BINARY_DIR}/contrib/gmsh
    PATH_SUFFIXES
    include include/gmsh
    DOC "Directory where GMSH header files are stored" )

  if ( GMSH_INCLUDE_PATH )
    set( FEELPP_HAS_GMSH_H 1 )
    FIND_PATH(GMSH_ADAPTMESH_INCLUDE_DIR
      Openfile.h Field.h
      PATHS ${GMSH_INCLUDE_PATH}
      DOC "Directory where GMSH header files are stored" )
    if ( GMSH_ADAPTMESH_INCLUDE_DIR )
      set( FEELPP_HAS_GMSH_H 1 )
    else ( GMSH_ADAPTMESH_INCLUDE_DIR )
      message(STATUS "Gmsh headers: some headers needed for meshadaptation are missing")
      message(STATUS "Check wiki pages for mesh adaptation to install properly gmsh")
    endif( GMSH_ADAPTMESH_INCLUDE_DIR )
  endif(GMSH_INCLUDE_PATH)

  FIND_LIBRARY(GMSH_LIBRARY NAMES Gmsh gmsh-2.5.1 gmsh1 gmsh
    PATHS
    $ENV{GMSH_DIR}
    ${CMAKE_BINARY_DIR}/contrib/gmsh
    ${CMAKE_SYSTEM_PREFIX_PATH}
    PATH_SUFFIXES
    lib lib/x86_64-linux-gnu/ )

  if( NOT GMSH_LIBRARY )
    if(APPLE)
      set( GMSHLIB libGmsh.dylib )
    else(APPLE)
      set( GMSHLIB libGmsh.so )
    endif(APPLE)
    FIND_PATH(GMSH_LIBRARY_PATH ${GMSHLIB}
      PATHS
      $ENV{GMSH_DIR}/lib
      ${CMAKE_BINARY_DIR}/contrib/gmsh/lib
      NO_DEFAULT_PATH)

    if(GMSH_LIBRARY_PATH)
        set(GMSH_LIBRARY "${GMSH_LIBRARY_PATH}/${GMSHLIB}" )
    else(GMSH_LIBRARY_PATH)
        set(GMSH_LIBRARY "GMSH_LIBRARY-NOTFOUND" )
    endif(GMSH_LIBRARY_PATH)
  endif(NOT GMSH_LIBRARY)

  FIND_LIBRARY(GL2PS_LIBRARY NAMES gl2ps
    PATH
    $ENV{GMSH_DIR}
    ${CMAKE_BINARY_DIR}/contrib/gmsh/lib
    ${CMAKE_SYSTEM_PREFIX_PATH}
    PATH_SUFFIXES
    lib  )

  IF ( FEELPP_ENABLE_OPENGL )
    FIND_LIBRARY(GL_LIBRARY NAMES GL
      PATH
      $ENV{GMSH_DIR}
      ${CMAKE_BINARY_DIR}/contrib/gmsh/
      PATH_SUFFIXES
      lib  )
  ENDIF()

  if(GMSH_INCLUDE_PATH)
      set(GMSH_INCLUDE_DIR ${GMSH_INCLUDE_PATH})
  endif(GMSH_INCLUDE_PATH)

  set(GMSH_LIBRARIES "")
  if(GMSH_LIBRARY)
    set(GMSH_LIBRARIES ${GMSH_LIBRARIES} ${GMSH_LIBRARY})
endif()

  FIND_PACKAGE_HANDLE_STANDARD_ARGS (GMSH DEFAULT_MSG
      GMSH_INCLUDE_DIR GMSH_LIBRARIES GMSH_EXECUTABLE
    )

  if ( GMSH_FOUND )
    set(FEELPP_HAS_GMSH_LIBRARY 1)
    MESSAGE( STATUS "GMSH found: header(${GMSH_INCLUDE_DIR}) lib(${GMSH_LIBRARY}) executable(${GMSH_EXECUTABLE})" )
    MESSAGE( STATUS "GL2PS found: lib(${GL2PS_LIBRARY})" )
    IF ( FEELPP_ENABLE_OPENGL )
      MESSAGE( STATUS "GL found: lib(${GL_LIBRARY})" )
    ENDIF()
  endif()

  mark_as_advanced( GMSH_INCLUDE_DIR )
  mark_as_advanced( GMSH_LIBRARY )
  mark_as_advanced( GL2PS_LIBRARY )
  IF ( FEELPP_ENABLE_OPENGL )
    mark_as_advanced( GL_LIBRARY )
  ENDIF()
  mark_as_advanced( GMSH_EXECUTABLE )

else(FEELPP_ENABLE_GMSH_LIBRARY)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS (GMSH DEFAULT_MSG GMSH_EXECUTABLE )

  if ( GMSH_FOUND )
    MESSAGE( STATUS "GMSH found: executable(${GMSH_EXECUTABLE})" )
  endif()
  mark_as_advanced( GMSH_EXECUTABLE )

endif(FEELPP_ENABLE_GMSH_LIBRARY)

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

find_program( GMSH_EXECUTABLE gmsh DOC "GMSH mesh generator" )


option(FEELPP_ENABLE_GMSH_LIBRARY "Enables Gmsh library in Feel++" ON )
if ( FEELPP_ENABLE_GMSH_LIBRARY )
  INCLUDE(CheckIncludeFileCXX)
  FIND_PATH(GMSH_INCLUDE_DIR
    Gmsh.h Context.h GModel.h
    PATHS ${CMAKE_SYSTEM_PREFIX_PATH} $ENV{GMSH_DIR}/include/gmsh
    PATH_SUFFIXES include include/gmsh
    DOC "Directory where GMSH header files are stored" )
  include_directories(${GMSH_INCLUDE_DIR})
  if ( GMSH_INCLUDE_DIR )
    FIND_PATH(GMSH_ADAPTMESH_INCLUDE_DIR
      Openfile.h Field.h
      PATHS ${GMSH_INCLUDE_DIR}
      DOC "Directory where GMSH header files are stored" )
      if ( GMSH_ADAPTMESH_INCLUDE_DIR )
	set( FEELPP_HAS_GMSH_H 1 )
      else ( GMSH_ADAPTMESH_INCLUDE_DIR )
	message(STATUS "Gmsh headers: some headers needed for meshadaptation are missing")
	message(STATUS "Check wiki pages for mesh adaptation to install properly gmsh")
      endif( GMSH_ADAPTMESH_INCLUDE_DIR )
  endif()
  #include(CheckIncludeFiles)
  #set(CMAKE_REQUIRED_INCLUDES "${GMSH_INCLUDE_DIR};${CMAKE_REQUIRED_INCLUDES}")
  #check_include_file(Gmsh.h FEELPP_HAS_GMSH_GMSH_H )
  ##check_include_file(Context.h FEELPP_HAS_GMSH_CONTEXT_H )
  #check_include_file(GModel.h FEELPP_HAS_GMSH_GMODEL_H )
  #if ( FEELPP_HAS_GMSH_GMODEL_H AND FEELPP_HAS_GMSH_CONTEXT_H and FEELPP_HAS_GMSH_GMSH_H )
  #  set( FEELPP_HAS_GMSH_H 1 )
  #endif()
  #message(STATUS "Gmsh headers : ${FEELPP_HAS_GMSH_H}, ${CMAKE_REQUIRED_INCLUDES}" )

  FIND_LIBRARY(GMSH_LIBRARY NAMES Gmsh gmsh-2.5.1 gmsh1 gmsh
    PATH
    ${CMAKE_SYSTEM_PREFIX_PATH}
    $ENV{GMSH_DIR}
    PATH_SUFFIXES
    lib  )
  if( NOT GMSH_LIBRARY )
    FIND_PATH(GMSH_LIBRARY_PATH
      libGmsh.so
      PATHS ${CMAKE_SYSTEM_PREFIX_PATH} $ENV{GMSH_DIR}/
      PATH_SUFFIXES lib )
    set(GMSH_LIBRARY "${GMSH_LIBRARY_PATH}/libGmsh.so" )
  endif()

  FIND_LIBRARY(GL2PS_LIBRARY NAMES gl2ps
    PATH
    ${CMAKE_SYSTEM_PREFIX_PATH}
    $ENV{GMSH_DIR}/lib
    PATH_SUFFIXES
    lib  )
  FIND_LIBRARY(GL_LIBRARY NAMES GL
    PATH
    ${CMAKE_SYSTEM_PREFIX_PATH}
    $ENV{GMSH_DIR}/lib
    PATH_SUFFIXES
    lib  )

  FIND_PACKAGE_HANDLE_STANDARD_ARGS (GMSH DEFAULT_MSG
    GMSH_INCLUDE_DIR GMSH_LIBRARY GMSH_EXECUTABLE
    )

  if ( GMSH_FOUND )
    set(FEELPP_HAS_GMSH_LIBRARY 1)
    MESSAGE( STATUS "GMSH found: header(${GMSH_INCLUDE_DIR}) lib(${GMSH_LIBRARY}) executable(${GMSH_EXECUTABLE})" )
    MESSAGE( STATUS "GL2PS found: lib(${GL2PS_LIBRARY})" )
    MESSAGE( STATUS "GL found: lib(${GL_LIBRARY})" )
  endif()

  mark_as_advanced( GMSH_INCLUDE_DIR )
  mark_as_advanced( GMSH_LIBRARY )
  mark_as_advanced( GL2PS_LIBRARY )
  mark_as_advanced( GL_LIBRARY )
  mark_as_advanced( GMSH_EXECUTABLE )

else(FEELPP_ENABLE_GMSH_LIBRARY)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS (GMSH DEFAULT_MSG GMSH_EXECUTABLE )

  if ( GMSH_FOUND )
    MESSAGE( STATUS "GMSH found: executable(${GMSH_EXECUTABLE})" )
  endif()
  mark_as_advanced( GMSH_EXECUTABLE )

endif(FEELPP_ENABLE_GMSH_LIBRARY)






# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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

function(_gmsh_get_version _out_major _out_minor _out_patch _gmsh_version_h)
	file(STRINGS ${_gmsh_version_h} _gmsh_vinfo REGEX "^#define[\t ]+GMSH_.*__VERSION.*")
	if (NOT _gmsh_vinfo)
		message(FATAL_ERROR "include file ${_gmsh_version_h} does not exist")
	endif()
	string(REGEX REPLACE "^.*GMSH_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_major} "${_gmsh_vinfo}")
	string(REGEX REPLACE "^.*GMSH_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_minor} "${_gmsh_vinfo}")
	string(REGEX REPLACE "^.*GMSH_PATCH_VERSION[ \t]+([0-9]+).*" "\\1" ${_out_patch} "${_gmsh_vinfo}")
	if (NOT ${_out_major} MATCHES "[0-9]+")
		message(FATAL_ERROR "failed to determine GMSH_MAJOR_VERSION, "
			            "expected a number, got ${${_out_major}}")
	endif()
	if (NOT ${_out_minor} MATCHES "[0-9]+")
		message(FATAL_ERROR "failed to determine GMSH_MINOR_VERSION, "
			            "expected a number, got ${${_out_minor}}")
	endif()
	if (NOT ${_out_patch} MATCHES "[0-9]+")
		message(FATAL_ERROR "failed to determine GMSH_PATCH_VERSION, "
			            "expected a number, got ${${_out_patch}}")
	endif()
	message(STATUS "found Gmsh [${_gmsh_version_h}], version ${${_out_major}}.${${_out_minor}}.${${_out_patch}}")
	set(${_out_major} ${${_out_major}} PARENT_SCOPE)
	set(${_out_minor} ${${_out_minor}} PARENT_SCOPE)
	set(${_out_patch} ${${_out_patch}} PARENT_SCOPE)
endfunction()


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
  DOC "GMSH mesh generator"
  NO_DEFAULT_PATH
  )
if(NOT GMSH_EXECUTABLE)
find_program( GMSH_EXECUTABLE
  NAMES gmsh
  PATHS
  $ENV{GMSH_DIR}/bin
  ${CMAKE_BINARY_DIR}/contrib/gmsh/bin
  PATH_SUFFIXES bin
  DOC "GMSH mesh generator" )
endif()

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
      set( FEELPP_HAS_GMSH_ADAPT_H 1 )
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

  OPTION( FEELPP_ENABLE_GL2PS "Enable the GL2PS library" ON )
  IF ( FEELPP_ENABLE_GL2PS )
    FIND_LIBRARY(GL2PS_LIBRARY NAMES gl2ps
      PATH
      $ENV{GMSH_DIR}
      ${CMAKE_BINARY_DIR}/contrib/gmsh/lib
      ${CMAKE_SYSTEM_PREFIX_PATH}
      PATH_SUFFIXES
      lib  )
  ENDIF( FEELPP_ENABLE_GL2PS )

  IF ( FEELPP_ENABLE_OPENGL )
    FIND_LIBRARY(GL_LIBRARY NAMES GL
      PATH
      $ENV{GMSH_DIR}
      ${CMAKE_BINARY_DIR}/contrib/gmsh/
      PATH_SUFFIXES
      lib  )
  ENDIF()

  # should detect instead if MED is supported by Gmsh
  OPTION( FEELPP_ENABLE_Med "Enable the MED library" OFF )
  OPTION( FEELPP_ENABLE_OpenCascade "Enable the OCC library" OFF )
  IF(GMSH_EXECUTABLE)
    execute_process(COMMAND ${GMSH_EXECUTABLE} -info
      OUTPUT_FILE "gmsh-info.log"
      ERROR_FILE "gmsh-config.err")
    execute_process(COMMAND grep "Build options" gmsh-config.err
      OUTPUT_VARIABLE GMSH_BUILDOPTS)
    #MESSAGE(STATUS "GMSH_BUILDOPTS=${GMSH_BUILDOPTS}")
    STRING(REPLACE "Build options    : " "" GMSH_EXTERNAL_LIBS_STRING ${GMSH_BUILDOPTS})
    string(REPLACE " " ";" GMSH_EXTERNAL_LIBS ${GMSH_EXTERNAL_LIBS_STRING})
    #MESSAGE(STATUS "GMSH_EXTERNAL_LIBS=${GMSH_EXTERNAL_LIBS}")

    set (EXTERNAL_LIBS Med OpenCascade)
    foreach(_l ${EXTERNAL_LIBS})
      list (FIND GMSH_EXTERNAL_LIBS ${_l} _index)
      if (${_index} GREATER -1)
        MESSAGE(STATUS "Gmsh has ${_l} support")
        set(FEELPP_ENABLE_${_l} ON)
      ENDIF()
    endforeach()
  ELSE(GMSH_EXECUTABLE)
    MESSAGE(STATUS "gmsh executable has not been detected. Support for Med and OpenCascade set to OFF")
  ENDIF(GMSH_EXECUTABLE)

  set(GMSH_EXTERNAL_LIBRARIES "")
  IF ( FEELPP_ENABLE_MED )
    find_library(MED_LIB med)
    if(MED_LIB)
      MESSAGE(WARNING "[Feelpp] Add support for MED library from Gmsh (requires HDF5)")
      list(APPEND GMSH_EXTERNAL_LIBRARIES ${MED_LIB})
    endif(MED_LIB)
    find_library(MEDC_LIB medC)
    if(MEDC_LIB)
      list(APPEND GMSH_EXTERNAL_LIBRARIES ${MEDC_LIB})
    endif(MEDC_LIB)
    if (MED_LIB AND MEDC_LIB)
      # check if HDF5 version is compatible
      IF (NOT HDF_VERSION_MAJOR_REF EQUAL 1 OR NOT HDF_VERSION_MINOR_REF EQUAL 8)
        MESSAGE(STATUS "[feelpp] HDF5 version is ${HDF_VERSION_REF}. Only 1.8.x versions are compatible with med support.")
      ELSE()
        ADD_DEFINITIONS( -DFEEPP_HAS_GMSH_HAS_MED )
      ENDIF()
    endif(MED_LIB AND MEDC_LIB)
  endif(FEELPP_ENABLE_MED)
  MESSAGE(STATUS "[Gmsh] GMSH_EXTERNAL_LIBRARIES ${GMSH_EXTERNAL_LIBRARIES}")
  
  if(GMSH_INCLUDE_PATH)
      set(GMSH_INCLUDE_DIR ${GMSH_INCLUDE_PATH})
  endif(GMSH_INCLUDE_PATH)

  set(GMSH_LIBRARIES "")
  if(GMSH_LIBRARY)
    set(GMSH_LIBRARIES ${GMSH_LIBRARIES} ${GMSH_LIBRARY} ${GMSH_EXTERNAL_LIBRARIES})
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

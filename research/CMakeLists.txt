###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-12
#
#  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)
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


project( feelpp-research )
feelpp_get_environment()

if ( (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) OR (FEELPP_COMPONENT STREQUAL "research" ))
  find_package(Feel++ COMPONENTS Toolboxes MOR PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[research] found Feel++")
  endif()
  if(NOT FEELPP_TOOLBOXES_FOUND)
    message(FATAL_ERROR "Feel++ Toolboxes were not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[research] found Feel++ toolboxes")
  endif()
  if(NOT FEELPP_MOR_FOUND)
    message(FATAL_ERROR "Feel++ MOR were not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[research] found Feel++ MOR")
  endif()

  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
  
  message(STATUS "[install] ${CMAKE_INSTALL_PREFIX}")
endif()

#
# list of sub-directories of curdir
MACRO(project_list_subdir result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

# Load documented project for Doxygen
# it assumes that a project is valid if :
#   there is README.{md,adoc} file in the directory
#   that the directory path include "research"
# for each detected project it creates an FEELPP_ENABLE_RESEARCH_${PROJECT} and FEELPP_ENABLE_RESEARCH_${PROJECT}_DOCUMENTATION
# option which are set to OFF by default
MESSAGE(STATUS "[feelpp] Feel++ projects :")
project_list_subdir(projects ${CMAKE_CURRENT_SOURCE_DIR})
foreach(project ${projects})
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${project}/README.md OR 
      EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${project}/README.adoc  )
    string(TOUPPER ${project} PROJECT)
    if ( ${PROJECT} STREQUAL "EYE2BRAIN" OR ${PROJECT} STREQUAL "BRAIN" OR ${PROJECT} STREQUAL "SWIMMER"  )
      option(FEELPP_ENABLE_RESEARCH_${PROJECT} "Enable/Disable research/${project} " ON)
    else()
      option(FEELPP_ENABLE_RESEARCH_${PROJECT} "Enable/Disable research/${project} " ON)
    endif()
    option(FEELPP_ENABLE_RESEARCH_${PROJECT}_DOCUMENTATION "Enable/Disable documentation in research/${project}" OFF)
		SET(FEELPP_ENABLED_PROJECTS "${FEELPP_ENABLED_PROJECTS} ${project}" )
		MESSAGE(STATUS "[feelpp] ${project}\t-- Compile : ${FEELPP_ENABLE_RESEARCH_${PROJECT}} ")
		MESSAGE(STATUS "        \t\t-- Doc :     ${FEELPP_ENABLE_RESEARCH_${PROJECT}_DOCUMENTATION}" )
  endif()
endforeach()

separate_arguments(FEELPP_ENABLED_PROJECTS)
FOREACH( project ${FEELPP_ENABLED_PROJECTS} )
	string(TOUPPER ${project} PROJECT)
    if ( ${project} STREQUAL "eye2brain" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_EYE2BRAIN ON )
    endif()
    if ( ${project} STREQUAL "hydromorpho" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_HYDROMORPHO OFF )
    endif()
    if ( ${project} STREQUAL "sivibirpp" AND  EXISTS ${project} )
      SET( FEELPP_ENABLE_RESEARCH_SIVIBIRPP OFF )
    endif()
	IF(${FEELPP_ENABLE_RESEARCH_${PROJECT}} )
		add_subdirectory( ${project} )
	endif()
endforeach()

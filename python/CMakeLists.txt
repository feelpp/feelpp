
###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-12
#
#  Copyright (C) 2009-2017 Universit� Joseph Fourier (Grenoble I)
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
cmake_minimum_required(VERSION 3.13.0)
if(DEFINED PROJECT_NAME)
  set( PYTHON_SUBPROJECT ON )
endif()

project(feelpp-python)
  
set(FEELPP_PYTHON_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(FEELPP_PYTHON_BINARY_DIR ${PROJECT_BINARY_DIR})

# with this it is possible to compile applications/MOR from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "python" )
  find_package(Feel++ COMPONENTS Toolboxes MOR PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "[python] Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[python] found Feel++ library")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(WARNING "[python] Feel++ was not installed with toolboxes component")
  else()
    message(STATUS "[python] found Feel++ toolboxes")
  endif()
    if (NOT FEELPP_MOR_FOUND)
    message(WARNING "[python] Feel++ was not installed with mor component")
  else()
    message(STATUS "[python] found Feel++ toolboxes")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  feelpp_get_environment()
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_TOOLBOXES_FOUND ON)
  set(FEELPP_MOR_FOUND ON)
endif()


if ( FEELPP_HAS_PYTHON )
  message(STATUS "[feelpp-python] python module path:${FEELPP_PYTHON_MODULE_PATH}")
endif()

#include(feelpp.module.pyfeelpp)

add_subdirectory(pyfeelpp)

# TOOLBOXES needs to be defined for pyfeelpp-toolboxes 
add_subdirectory( pyfeelpp-toolboxes )

add_subdirectory(pyfeelpp-mor)
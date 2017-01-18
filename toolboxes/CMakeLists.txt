###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 05 Jan 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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

if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} GREATER 3 )
  message(STATUS "[feelpp] use OLD policy CMP0064" )
  cmake_policy(SET CMP0064 OLD)
endif()

# with this it is possible to compile applications/models from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feel/cmake/modules /usr/share/feel/cmake/modules /usr/local/share/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

# if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
#   cmake_minimum_required(VERSION 2.8)
#   FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH $ENV{FEELPP_DIR}/share/feel/cmake/modules /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
#   if ( FEELPP_CMAKE_MODULES )
#     set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
#   else()
#     message(FATAL_ERROR "Feel++ does not seem to have been installed on this platform")
#   endif()
#   #OPTION(FEELPP_ENABLE_SYSTEM_EIGEN3 "enable system eigen3 support" ON)
#   Find_Package(Feel++)
#   include(feelpp.macros)
# endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
if(1)#FEELPP_HAS_GMSH)
  add_subdirectory( feel/feelmodels )
  add_subdirectory( thermodyn )
  add_subdirectory( fluid )
  add_subdirectory( solid )
  add_subdirectory( fsi )
  add_subdirectory( hdg )
  #add_subdirectory( advection )
  add_subdirectory( levelset )
  add_subdirectory( multifluid )
  add_subdirectory( thermoelectric )
endif()

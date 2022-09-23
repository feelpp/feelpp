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
  set( MOR_SUBPROJECT ON )
endif()

project(mor)
  

set(FEELPP_MOR_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(FEELPP_MOR_BINARY_DIR ${PROJECT_BINARY_DIR})

# with this it is possible to compile applications/MOR from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "mor" )
  find_package(Feel++ COMPONENTS Toolboxes PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(WARNING "[mor] Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[mor] found Feel++ library")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(WARNING "[mor] Feel++ was not installed with toolboxes component")
  else()
    message(STATUS "[mor] found Feel++ toolboxes")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  feelpp_get_environment()
  # useful for pyfeelpp-toolboxes
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
#  include_directories(${FEELPP_MOR_SOURCE_DIR})
#  include_directories(${FEELPP_MOR_BINARY_DIR})  
  set(FEELPP_FOUND ON)  
  set(FEELPP_TOOLBOXES_FOUND ON)
endif()



OPTION(FEELPP_ENABLE_MOR_ALL "enable all MOR testcase" OFF)
OPTION(FEELPP_ENABLE_MOR_EXAMPLES "enable MOR examples" ON)
OPTION(FEELPP_ENABLE_MOR_MORDICUS "enable mordicus" OFF)


set(FEELPP_MOR_LIBRARIES)

INCLUDE(CTest)
ENABLE_TESTING()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
include(feelpp.macros.mor)
add_subdirectory(cmake)

if ( FEELPP_FOUND AND FEELPP_TOOLBOXES_FOUND )
  add_subdirectory(mor/feel/feelmor)
  add_subdirectory(apps)




if ( FEELPP_ENABLE_MOR_MORDICUS )
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mordicus/src/cpp )
    set(MORDICUS_ENABLE_CEMOSIS ON)
    add_subdirectory(mordicus)
    message(STATUS "[feelpp-mor] Mordicus enabled")
  else()
    set(FEELPP_ENABLE_MOR_MORDICUS OFF)
    message(WARNING "[feelpp-mor] Mordicus enabled but not present")
  endif()
ENDIF()


OPTION(FEELPP_ENABLE_TESTS "enable Feel++ tests" ON)
if (FEELPP_ENABLE_TESTS  )
  add_subdirectory(tests)
endif()

OPTION(FEELPP_ENABLE_EXAMPLES "enable Feel++ examples" OFF)
if ( FEELPP_ENABLE_MOR_EXAMPLES )
  add_subdirectory( examples )
endif()


include(${FEELPP_MOR_SOURCE_DIR}/cmake/feelpp.mor.cmake)
install(EXPORT feelpp-mor-export-targets
  FILE feelpp-mor-export-targets.cmake
  NAMESPACE Feelpp::
  DESTINATION share/feelpp/mor/cmake/
  )

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp-mor.package.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp.version.cmake)
endif()

endif(FEELPP_FOUND AND FEELPP_TOOLBOXES_FOUND)
###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 19 Jun 2017
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
#
# setup pyfeelpp-mor
cmake_minimum_required(VERSION 3.13)
project(pyfeelpp-mor)

if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "mor/pyfeelpp-mor" )
  find_package(Feel++ COMPONENTS Toolboxes MOR PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "[pyfeelpp-mor] Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[pyfeelpp-mor] found Feel++ library")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(WARNING "[pyfeelpp-mor] Feel++ was not installed with toolboxes component")
  else()
    message(STATUS "[pyfeelpp-mor] found Feel++ toolboxes")
  endif()
  if (NOT FEELPP_MOR_FOUND)
    message(WARNING "[pyfeelpp-mor] Feel++ was not installed with mor component")
  else()
    message(STATUS "[pyfeelpp-mor] found Feel++ mor")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
endif()


OPTION(FEELPP_ENABLE_PYFEELPP_MOR "Enable PyFeel++ MOR" ON  )
if (FEELPP_ENABLE_PYFEELPP_MOR)
  set(PYFILES generate_basis.py )
  install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/mor)

  message(STATUS "Enabling pyfeelpp mor")
  add_subdirectory(feelpp)
  add_subdirectory(tests)
endif()

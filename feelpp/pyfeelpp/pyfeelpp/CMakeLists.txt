###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 23 Jun 2017
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
cmake_minimum_required(VERSION 3.3)
project(pyfeelpp)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

enable_testing()

if ( NOT TARGET check )
  add_custom_target(check)
endif()

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  list(APPEND CMAKE_MODULE_PATH "$ENV{FEELPP_DIR}/share/cmake/pybind11")
  list(APPEND CMAKE_MODULE_PATH "/usr/local/share/cmake/pybind11")
  include(pybind11Tools)
endif()

pybind11_add_module(_pyfeelpp python.cpp)
target_include_directories(_pyfeelpp PRIVATE ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS})
target_link_libraries( _pyfeelpp PRIVATE Feelpp::feelpp )

#add_subdirectory(core)

message(STATUS "[pybind11] ${PYBIND11_INCLUDE_DIRS}")
pybind11_add_module(test_nested  test_nested.cpp)
target_include_directories(test_nested PRIVATE ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS})
target_link_libraries( test_nested PRIVATE Feelpp::feelpp pybind11::pybind11 )




#add_subdirectory(filters)
#add_subdirectory(toolboxes)

#configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)





#add_subdirectory(core)
#add_subdirectory(crb)

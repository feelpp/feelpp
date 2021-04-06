###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 06 Aug 2018
#
#  Copyright (C) 2018 Feel++ Consortium
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
project(toolboxes-thermoelectric)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ COMPONENTS Toolboxes PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()


pybind11_add_module(_thermoelectric  thermoelectric.cpp)
target_include_directories(_thermoelectric PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR} ${FEELPP_DIR}/include/feelpp/toolboxes)
find_library( FEELPP_THERMOELECTRIC feelpp_toolbox_thermoelectric_lib HINTS ${FEELPP_DIR}/lib)
target_link_libraries( _thermoelectric PRIVATE ${FEELPP_THERMOELECTRIC}  ${FEELPP_LIBRARIES} )

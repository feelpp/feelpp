###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-28
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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
project(testsuite)
set(PROJECT_SHORTNAME "test")

if ( ${CMAKE_SOURCE_DIR}/ STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
else()
  feelpp_get_environment()
endif()

enable_testing()

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_custom_target( testsuite )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#add_subdirectory(contrib)

option(FEELPP_DISABLE_COMPILE_FAILURE "disable compilation failure in testsuite" OFF)
if (FEELPP_DISABLE_COMPILE_FAILURE )
  add_definitions(-DFEELPP_DISABLE_COMPILE_FAILURE=1)
endif()
message(STATUS "[feelpp/testsuite] disable compile failure: ${FEELPP_DISABLE_COMPILE_FAILURE}")



foreach( testdir core math alg mesh poly filters material models interpolation integration discr crb ts vf leaks opt ls pde fit fmi)
  add_custom_target( test${testdir} )
  add_subdirectory( feel${testdir} )

endforeach()

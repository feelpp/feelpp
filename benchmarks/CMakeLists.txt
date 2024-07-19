###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 07 Apr 2019
#
#  Copyright (C) 2019 Feel++ Consortium
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


cmake_minimum_required(VERSION 3.13)

project(benchmarks)
set(PROJECT_SHORTNAME "bench")

# with this it is possible to compile applications/models from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR ) OR ( FEELPP_COMPONENT STREQUAL "benchmarks" ) )
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
else()
  feelpp_get_environment()
  # usefull for pyfeelpp-toolboxes
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
endif()

option(BENCHMARK_ENABLE_TESTING  "Enable building the unit tests which depend on gtest" OFF)
set( BENCHMARK_ENABLE_TESTING OFF )
feelppGitSubmodulePrepare( benchmark )
include_directories(${CMAKE_SOURCE_DIR}/benchmarks/benchmark/include)
add_subdirectory(benchmark)
add_subdirectory(feelpp)
add_subdirectory(compilation)
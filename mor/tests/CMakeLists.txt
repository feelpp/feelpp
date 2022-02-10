###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2022-02-05
#
#  Copyright (C) 2022 Université de Strasbourg
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
project(testsuite-mor)
set(PROJECT_SHORTNAME "test")

enable_testing()

if ( NOT TARGET testsuite )
    add_custom_target( testsuite )
endif()

option(FEELPP_DISABLE_COMPILE_FAILURE "disable compilation failure in testsuite" OFF)
if (FEELPP_DISABLE_COMPILE_FAILURE )
  add_definitions(-DFEELPP_DISABLE_COMPILE_FAILURE=1)
endif()
message(STATUS "[feelpp/testsuite] disable compile failure: ${FEELPP_DISABLE_COMPILE_FAILURE}")



foreach( testdir mor)
  add_custom_target( test${testdir} )
  add_subdirectory( feel${testdir} )
endforeach()

###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 10 Jul 2017
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
# PyBind11
#
OPTION( FEELPP_ENABLE_PYBIND11 "Enable PYBIND11" ON )

if ( FEELPP_ENABLE_PYBIND11 )
  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/pybind11 )
    if ( GIT_FOUND  AND EXISTS ${CMAKE_SOURCE_DIR}/.git )
      execute_process(
        COMMAND git submodule update --init --recursive contrib/pybind11
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_FILE git.pybind11.log
        ERROR_FILE git.pybind11.log
        RESULT_VARIABLE ERROR_CODE
        )
      if(ERROR_CODE EQUAL "0")
        MESSAGE(STATUS "[feelpp] Git submodule contrib/pybind11 updated.")
      else()
        MESSAGE(WARNING "Git submodule contrib/pybind11 failed to be updated. Possible cause: No internet access, firewalls ...")
      endif()
    else()
      if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/pybind11/ )
        message( WARNING "Please make sure that git submodule contrib/pybind11 is available")
        message( WARNING "  run `git submodule update --init --recursive contrib/pybind11`")
      endif()
    endif()
    
  endif()


  add_subdirectory(contrib/pybind11)
  include_directories(${PYBIND11_INCLUDE_DIR})
  SET(FEELPP_HAS_PYBIND11 1)
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PyBind11" )
  ADD_DEFINITIONS( -DFEELPP_HAS_PYBIND11 )
endif()

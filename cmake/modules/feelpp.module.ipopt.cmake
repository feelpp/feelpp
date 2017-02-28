###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-19
#
#  Copyright (C) 2014-2015 Feel++ Consortium
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
OPTION( FEELPP_ENABLE_IPOPT "Enable IPOPT (Interior Point OPTimizer Library)" OFF )

if ( FEELPP_ENABLE_IPOPT )

  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/ipopt )

    if ( GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git )

      execute_process(
        COMMAND git submodule update --init --recursive contrib/ipopt
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_FILE ${CMAKE_BUILD_DIR}/git.ipopt.log
        ERROR_FILE ${CMAKE_BUILD_DIR}/git.ipopt.log
        RESULT_VARIABLE ERROR_CODE
        )
      if(ERROR_CODE EQUAL "0")
        MESSAGE(STATUS "Git submodule contrib/ipopt updated.")
      else()
        MESSAGE(WARNING "Git submodule contrib/ipopt failed to be updated. Possible cause: No internet access, firewalls ...")
      endif()
    else()
      if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/ )
        message( WARNING "Please make sure that git submodule contrib/ipopt is available")
        message( WARNING "  run `git submodule update --init --recursive contrib/ipopt`")
      endif()
    endif()

    if ( EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/CMakeLists.txt )

      message(STATUS "[feelpp] use contrib/ipopt : ${CMAKE_SOURCE_DIR}/contrib/ipopt/")
      SET(FEELPP_HAS_IPOPT 1)
      #ADD_DEFINITIONS( -DFEELPP_HAS_IPOPT )
      #ADD_DEFINITIONS( -fPIC )

      SET(FEELPP_LIBRARIES feelpp_ipopt ${FEELPP_LIBRARIES})
      SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} Ipopt/Contrib" )
      #add_subdirectory(contrib/ipopt)

      #SET(IPOPT_INCLUDE_DIR
      #          ${FEELPP_SOURCE_DIR}/contrib/ipopt/Ipopt/src/Interfaces
      #)

      # Compile/copy header in cmake binary dirs.
      include_directories(${CMAKE_BINARY_DIR}/contrib/ipopt/include/)
    endif()
  endif()
endif()


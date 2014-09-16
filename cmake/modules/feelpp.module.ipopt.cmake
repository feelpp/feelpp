###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-19
#
#  Copyright (C) 2014 Feel++ Consortium
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

  if ( GIT_FOUND )

    execute_process(
      COMMAND git submodule update --init --recursive contrib/ipopt
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      )
    MESSAGE(STATUS "Git submodule contrib/ipopt updated.")
  else()
    if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/ )
      message( FATAL_ERROR "Please make sure that git submodule contrib/ipopt is available")
      message( FATAL_ERROR "  run `git submodule update --init --recursive contrib/ipopt`")
    endif()
  endif()

  # if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/api/ipopt.hpp )

  #   execute_process(
  #     COMMAND  touch  swig/ipopt.scm.in
  #     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/ipopt/
  #     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-touch"
  #     #ERROR_FILE "ipopt-touch-errors"
  #     )
  #   execute_process(
  #     COMMAND autoreconf --verbose --install --force
  #     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/ipopt/
  #     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-autoreconf"
  #     #ERROR_FILE "ipopt-autoreconf-errors"
  #     )
  #   execute_process(
  #     COMMAND autoreconf --verbose --install --force
  #     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/ipopt/
  #     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-autoreconf"
  #     #ERROR_FILE "ipopt-autoreconf-errors"
  #     )

  #   if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/configure )
  #     message(FATAL_ERROR "configure not available")
  #   endif()

  #   # ensure that build dir is created
  #   file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt)

  #   if ( FEELPP_ENABLE_BUILD_STATIC )
  #     execute_process(
  #       COMMAND ${FEELPP_HOME_DIR}/contrib/ipopt/configure --enable-maintainer-mode --with-cxx=yes CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
  #       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
  #       #OUTPUT_QUIET
  #       OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-configure"
  #       #ERROR_FILE "ipopt-configure-errors"
  #       )
  #   else()
  #     execute_process(
  #       COMMAND ${FEELPP_HOME_DIR}/contrib/ipopt/configure --enable-maintainer-mode --with-cxx=yes --enable-shared CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
  #       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
  #       #      OUTPUT_QUIET
  #       OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-configure"
  #       #ERROR_FILE "ipopt-configure-errors"
  #       )
  #   endif()
  #   execute_process(
  #     COMMAND make ipopt.hpp
  #     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt/api
  #     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-ipopthpp" )

  #   # delete all Makefiles before Cmake generate its own
  #   execute_process(
  #     COMMAND make distclean
  #     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
  #     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-distclean"
  #     )
  # endif()
  #if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/api/ipopt.hpp )
  #  message(FATAL_ERROR "NLOpt: ipopt.hpp was not generated")
  #else()
  #  message(STATUS "NLOpt: ipopt.hpp is generated")
  #endif()
  #include_directories(${FEELPP_SOURCE_DIR}/contrib/ipopt/api)
  #add_subdirectory(contrib/ipopt)
  #SET(FEELPP_LIBRARIES feelpp_ipopt ${FEELPP_LIBRARIES} )
  #SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} NLOpt" )
  #SET(FEELPP_HAS_IPOPT 1)

endif()

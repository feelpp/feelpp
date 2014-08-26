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
OPTION( FEELPP_ENABLE_NLOPT "Enable NLOPT (NonLinear Optimisation Library)" ON )

if ( FEELPP_ENABLE_NLOPT )

  if ( GIT_FOUND )

    execute_process(
      COMMAND git submodule update --init --recursive contrib/nlopt
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      )
    MESSAGE(STATUS "Git submodule contrib/nlopt updated.")
  else()
    if ( NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/ )
      message( FATAL_ERROR "Please make sure that git submodule contrib/nlopt is available")
      message( FATAL_ERROR "  run `git submodule update --init --recursive contrib/nlopt`")
    endif()
  endif()

  if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/api/nlopt.hpp )

    execute_process(
      COMMAND  touch  swig/nlopt.scm.in
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/nlopt/
      OUTPUT_FILE "nlopt-touch"
      #ERROR_FILE "nlopt-touch-errors"
      )
    execute_process(
      COMMAND autoreconf --verbose --install --force
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/nlopt/
      OUTPUT_FILE "nlopt-autoreconf"
      #ERROR_FILE "nlopt-autoreconf-errors"
      )
    execute_process(
      COMMAND autoreconf --verbose --install --force
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/nlopt/
      OUTPUT_FILE "nlopt-autoreconf"
      #ERROR_FILE "nlopt-autoreconf-errors"
      )

    if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/configure )
      message(FATAL_ERROR "configure not available")
    endif()

    # ensure that build dir is created
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt)

    if ( FEELPP_ENABLE_BUILD_STATIC )
      execute_process(
        COMMAND ${FEELPP_HOME_DIR}/contrib/nlopt/configure --enable-maintainer-mode --with-cxx=yes CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
        #OUTPUT_QUIET
        OUTPUT_FILE "nlopt-configure"
        #ERROR_FILE "nlopt-configure-errors"
        )
    else()
      execute_process(
        COMMAND ${FEELPP_HOME_DIR}/contrib/nlopt/configure --enable-maintainer-mode --with-cxx=yes --enable-shared CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
        #      OUTPUT_QUIET
        OUTPUT_FILE "nlopt-configure"
        #ERROR_FILE "nlopt-configure-errors"
        )
    endif()
    execute_process(
      COMMAND make nlopt.hpp
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt/api
      OUTPUT_FILE "nlopt-nlopthpp" )

    # delete all Makefiles before Cmake generate its own
    execute_process(
      COMMAND make distclean
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
      OUTPUT_FILE "nlopt-distclean"
      )
  endif()
  if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/api/nlopt.hpp )
    message(FATAL_ERROR "NLOpt: nlopt.hpp was not generated")
  else()
    message(STATUS "NLOpt: nlopt.hpp is generated")
  endif()
  include_directories(${FEELPP_SOURCE_DIR}/contrib/nlopt/api)
  add_subdirectory(contrib/nlopt)
  SET(FEELPP_LIBRARIES feelpp_nlopt ${FEELPP_LIBRARIES} )
  SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} NLOpt" )
  SET(FEELPP_HAS_NLOPT 1)

endif()

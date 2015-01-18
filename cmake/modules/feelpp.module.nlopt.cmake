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

  if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/nlopt )

    if ( GIT_FOUND )

      execute_process(
        COMMAND git submodule update --init --recursive contrib/nlopt
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_FILE git.nlopt.log
        ERROR_FILE git.nlopt.log
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
        OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-touch"
        #ERROR_FILE "nlopt-touch-errors"
        )
      execute_process(
        COMMAND autoreconf --verbose --install --force
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/nlopt/
        OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-autoreconf"
        #ERROR_FILE "nlopt-autoreconf-errors"
        )
      execute_process(
        COMMAND autoreconf --verbose --install --force
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/contrib/nlopt/
        OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-autoreconf"
        #ERROR_FILE "nlopt-autoreconf-errors"
        )

      if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/configure )
          message(FATAL_ERROR "configure not available (Possible cause: autotools might be missing from your system)")
      endif()

      # ensure that build dir is created
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt)

      if ( FEELPP_ENABLE_BUILD_STATIC )
        execute_process(
          COMMAND ${FEELPP_HOME_DIR}/contrib/nlopt/configure --enable-maintainer-mode --with-cxx=yes CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
          #OUTPUT_QUIET
          OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-configure"
          #ERROR_FILE "nlopt-configure-errors"
          )
      else()
        execute_process(
          COMMAND ${FEELPP_HOME_DIR}/contrib/nlopt/configure --enable-maintainer-mode --with-cxx=yes --enable-shared CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_CXX_COMPILER}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
          #      OUTPUT_QUIET
          OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-configure"
          #ERROR_FILE "nlopt-configure-errors"
          )
      endif()
      execute_process(
        COMMAND make nlopt.hpp
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt/api
        OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-nlopthpp" )

      # delete all Makefiles before Cmake generate its own
      execute_process(
        COMMAND make distclean
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
        OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-distclean"
        )
    endif()
    if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/nlopt/api/nlopt.hpp )
      message(FATAL_ERROR "NLOpt: nlopt.hpp was not generated")
    else()
      message(STATUS "NLOpt: nlopt.hpp is generated")
    endif()
    set(NLOPT_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/nlopt/api)
    add_subdirectory(contrib/nlopt)
    include_directories(${NLOPT_INCLUDE_DIR})
    SET(FEELPP_LIBRARIES feelpp_nlopt ${FEELPP_LIBRARIES} )
    SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} NLOpt" )
    SET(FEELPP_HAS_NLOPT 1)

  else( EXISTS ${CMAKE_SOURCE_DIR}/contrib/nlopt )

    FIND_PATH(NLOPT_INCLUDE_DIR nlopt.hpp
      $ENV{FEELPP_DIR}/include/feel/nlopt
      NO_DEFAULT_PATH)

    FIND_PATH(NLOPT_INCLUDE_DIR nlopt.hpp
      $ENV{FEELPP_DIR}/include/feel
      /usr/include/feel/nlopt
      /usr/local/include/feel/nlopt
      /opt/local/include/feel/nlopt
      NO_DEFAULT_PATH)
    message(STATUS "NLopt/system: ${NLOPT_INCLUDE_DIR}")
    include_directories(${NLOPT_INCLUDE_DIR})
    FIND_LIBRARY(NLOPT_LIBRARY  NAMES feelpp_nlopt  PATHS   $ENV{FEELPP_DIR}/lib  NO_DEFAULT_PATH)
    FIND_LIBRARY(NLOPT_LIBRARY  NAMES feelpp_nlopt    )

    string(REPLACE "include" "" NLOPT_DIR ${NLOPT_INCLUDE_DIR} )

    set(NLOPT_LIBRARIES ${NLOPT_LIBRARY})
    message(STATUS "[feelpp] loading nlopt from includes: ${NLOPT_INCLUDE_DIR} Libraries: ${NLOPT_LIBRARIES} Dir: ${NLOPT_DIR}" )


    # handle the QUIETLY and REQUIRED arguments and set NLOPT_FOUND to TRUE if
    # all listed variables are TRUE
    include (FindPackageHandleStandardArgs)
    find_package_handle_standard_args (NLOPT DEFAULT_MSG NLOPT_INCLUDE_DIR NLOPT_LIBRARIES NLOPT_DIR )

    mark_as_advanced (NLOPT_INCLUDE_DIR NLOPT_LIBRARIES NLOPT_DIR NLOPT_LIBRARY)

  endif(  EXISTS ${CMAKE_SOURCE_DIR}/contrib/nlopt  )

endif()

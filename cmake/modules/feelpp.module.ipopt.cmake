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
OPTION( FEELPP_ENABLE_IPOPT "Enable IPOPT (Interior Point OPTimizer)" ON )

if ( FEELPP_ENABLE_IPOPT )
    find_package( BLAS REQUIRED )
    if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/ipopt )
        if ( GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git )
            execute_process(
                COMMAND git submodule update --init --recursive contrib/ipopt
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_FILE git.ipopt.log
                ERROR_FILE git.ipopt.log
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
                message( WARNING "  run `git submodule update --init --recursive contrib/nlopt`")
            endif()
        endif()

        if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/Ipopt/src/ipopt.hpp )

            if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/configure )
                message(FATAL_ERROR "configure not available (Possible cause: autotools might be missing from your system)")
            endif()

            # ensure that build dir is created
            file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt)


            message(STATUS "IPOpt: configuring...")
            if ( FEELPP_ENABLE_BUILD_STATIC )
                execute_process(
                    COMMAND ${FEELPP_HOME_DIR}/contrib/ipopt/configure
                                --with-cxx=yes
                                --with-blas="-L/usr/lib/atlas-base/ -lf77blas -latlas -lcblas"
                                CXX=${CMAKE_CXX_COMPILER}
                                CC=${CMAKE_CXX_COMPILER}
                    #--enable-maintainer-mode
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
                    #OUTPUT_QUIET
                    OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-configure"
                    #ERROR_FILE "nlopt-configure-errors"
                    )
            else()
                execute_process(
                    COMMAND ${FEELPP_HOME_DIR}/contrib/ipopt/configure
                                --with-cxx=yes
                                --with-blas="-L/usr/lib/atlas-base/ -lf77blas -latlas -lcblas"
                                --enable-shared
                                CXX=${CMAKE_CXX_COMPILER}
                                CC=${CMAKE_CXX_COMPILER}
                    #--enable-maintainer-mode
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
                    #      OUTPUT_QUIET
                    OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-configure"
                    #ERROR_FILE "nlopt-configure-errors"
                    )
            endif()

            message(STATUS "IPOpt: compiling...")
            execute_process(
                COMMAND make
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/ipopt
                OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/ipopt-ipopthpp" )

            # delete all Makefiles before Cmake generate its own
            execute_process(
                COMMAND make distclean
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/contrib/nlopt
                OUTPUT_FILE "${CMAKE_SOURCE_DIR}/contrib/nlopt-distclean"
                )
        endif()

        if (NOT EXISTS ${FEELPP_SOURCE_DIR}/contrib/ipopt/Ipopt/src/ipopt.hpp )
            message(FATAL_ERROR "IPOpt: ipopt.hpp was not generated")
        else()
            message(STATUS "IPOpt: ipopt.hpp is generated")
        endif()

        set(IPOPT_INCLUDE_DIR ${FEELPP_SOURCE_DIR}/contrib/ipopt/Ipopt/src/)
        #add_subdirectory(contrib/nlopt)
        include_directories(${IPOPT_INCLUDE_DIR})
        SET(FEELPP_LIBRARIES feelpp_ipopt ${FEELPP_LIBRARIES} )
        SET(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} IPOpt" )
        SET(FEELPP_HAS_IPOPT 1)

    else( EXISTS ${CMAKE_SOURCE_DIR}/contrib/ipopt )

        FIND_PATH(IPOPT_INCLUDE_DIR ipopt.hpp
            $ENV{FEELPP_DIR}/include/feel/ipopt
            NO_DEFAULT_PATH)

        FIND_PATH(IPOPT_INCLUDE_DIR ipopt.hpp
            $ENV{FEELPP_DIR}/include/feel
            /usr/include/feel/ipopt
            /usr/local/include/feel/ipopt
            /opt/local/include/feel/ipopt
            NO_DEFAULT_PATH)
        message(STATUS "IPopt/system: ${IPOPT_INCLUDE_DIR}")
        include_directories(${IPOPT_INCLUDE_DIR})
        FIND_LIBRARY(IPOPT_LIBRARY  NAMES feelpp_ipopt  PATHS   $ENV{FEELPP_DIR}/lib  NO_DEFAULT_PATH)
        FIND_LIBRARY(IPOPT_LIBRARY  NAMES feelpp_ipopt    )

        string(REPLACE "include" "" IPOPT_DIR ${IPOPT_INCLUDE_DIR} )

        set(IPOPT_LIBRARIES ${IPOPT_LIBRARY})
        message(STATUS "[feelpp] loading ipopt from includes: ${IPOPT_INCLUDE_DIR} Libraries: ${IPOPT_LIBRARIES} Dir: ${IPOPT_DIR}" )

        if( NOT IPOPT_INCLUDE_DIR OR NOT IPOPT_LIBRARIES )
            message(FATAL_ERROR "IPopt was not found on your system. Either install it or set FEELPP_ENABLE_IPOPT to OFF.")
        endif()

        # handle the QUIETLY and REQUIRED arguments and set IPOPT_FOUND to TRUE if
        # all listed variables are TRUE
        include (FindPackageHandleStandardArgs)
        find_package_handle_standard_args (IPOPT DEFAULT_MSG IPOPT_INCLUDE_DIR IPOPT_LIBRARIES IPOPT_DIR )

        mark_as_advanced (IPOPT_INCLUDE_DIR IPOPT_LIBRARIES IPOPT_DIR IPOPT_LIBRARY)

    endif(  EXISTS ${CMAKE_SOURCE_DIR}/contrib/ipopt  )

endif()

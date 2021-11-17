###  feelpp.directive.cmake; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-04
#
#  Copyright (C) 2013-2015 Feel++ Consortium
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

set(FEELPP_SCHEDULER "None")
if ( FEELPP_ENABLE_SCHED_LOADLEVELER )
   set( FEELPP_SCHEDULER "LoadLeveler" )
elseif ( FEELPP_ENABLE_SCHED_SLURM )
     set( FEELPP_SCHEDULER "Slurm"  )
elseif ( FEELPP_ENABLE_SCHED_OAR )
     set( FEELPP_SCHEDULER "oar"  )
elseif ( FEELPP_ENABLE_SCHED_CCC )
     set( FEELPP_SCHEDULER "CCC"  )
endif()

################################################################################
# Display a post-configuration installation directives list
################################################################################

function(list_options varTarget varOut )
  set(res "")
  get_property(CD TARGET ${varTarget} PROPERTY INTERFACE_COMPILE_DEFINITIONS)
  message(STATUS "cd: ${CD} ")
  foreach( opts IN LISTS CD )
    string( REGEX MATCH "FEELPP_HAS_([a-zA-Z0-9]+)$" OPT ${opts} )
    #message( STATUS "match: ${CMAKE_MATCH_1}" )
    #message( STATUS "opt: ${OPT}" )
    if ( OPT )
      #message( STATUS "match: ${CMAKE_MATCH_1}" )
      list(APPEND res "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  list(REVERSE res )
  set(${varOut} ${res} PARENT_SCOPE )
endfunction()

function(feelpp_print_directive)
  list_options(Feelpp::feelpp_contrib FEELPP_CONTRIB_ENABLED_OPTIONS)
  list_options(Feelpp::feelpp FEELPP_ENABLED_OPTIONS)

  get_property(FEELPP_STD_CPP TARGET Feelpp::feelpp PROPERTY FEELPP_STD_CPP)

  MESSAGE(STATUS "================================================================================")
  MESSAGE(STATUS "     FEELPP_VERSION_MAJOR : ${FEELPP_VERSION_MAJOR}")
  MESSAGE(STATUS "     FEELPP_VERSION_MINOR : ${FEELPP_VERSION_MINOR}")
  MESSAGE(STATUS "     FEELPP_VERSION_MICRO : ${FEELPP_VERSION_MICRO}")
  MESSAGE(STATUS "FEELPP_VERSION_PRERELEASE : ${FEELPP_VERSION_PRERELEASE}")
  MESSAGE(STATUS "  FEELPP_VERSION_METADATA : ${FEELPP_VERSION_METADATA}")
  MESSAGE(STATUS "    FEELPP_VERSION_STRING : ${FEELPP_VERSION_STRING}")
  MESSAGE(STATUS "           FEELPP_REVISON : ${FEELPP_REVISION}")
  MESSAGE(STATUS "           FEELPP_BUILDID : ${FEELPP_BUILDID}")
  MESSAGE(STATUS "           FEELPP_STD_CPP : c++${FEELPP_STD_CPP}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Feel++ Modules :")
  MESSAGE(STATUS "     QuickStart: ${FEELPP_ENABLE_QUICKSTART}")
  MESSAGE(STATUS "  Documentation: ${FEELPP_ENABLE_DOCUMENTATION}")
  MESSAGE(STATUS "        Doxygen: ${FEELPP_ENABLE_DOXYGEN}")
  MESSAGE(STATUS "      Testsuite: ${FEELPP_ENABLE_TESTS}")
  MESSAGE(STATUS "   Applications: ${FEELPP_ENABLE_APPLICATIONS}")
  MESSAGE(STATUS "     Benchmarks: ${FEELPP_ENABLE_BENCHMARKS}")
  MESSAGE(STATUS "       Research: ${FEELPP_ENABLE_RESEARCH}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "Feel++ Configuration:")
  MESSAGE(STATUS "                CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
  MESSAGE(STATUS "              CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
  MESSAGE(STATUS "                 CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
  MESSAGE(STATUS "         CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
  MESSAGE(STATUS "  CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  MESSAGE(STATUS "           CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
  MESSAGE(STATUS "       CMAKE_CXX_FLAGS_DEBUGFULL: ${CMAKE_CXX_FLAGS_DEBUGFULL}")
  MESSAGE(STATUS "            CMAKE_CXX_FLAGS_ASAN: ${CMAKE_CXX_FLAGS_ASAN}")
  MESSAGE(STATUS "        CMAKE_CXX_FLAGS_COVERAGE: ${CMAKE_CXX_FLAGS_COVERAGE}")
  MESSAGE(STATUS "              FEELPP_INCLUDE_DIR: ${FEELPP_INCLUDE_DIR}")
  MESSAGE(STATUS "                 FEELPP_DATA_DIR: ${FEELPP_DATA_DIR}")
  MESSAGE(STATUS "            CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "             FEELPP_MACHINE_NAME: ${FEELPP_MACHINE_NAME}")
  MESSAGE(STATUS "               FEELPP_ENABLE_GIT: ${FEELPP_ENABLE_GIT}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "          FEELPP_ENABLE_MPI_MODE: ${FEELPP_ENABLE_MPI_MODE} Nprocs: ${N} Nprocs(build): ${N2}")
  MESSAGE(STATUS "               FEELPP_ENABLE_PCH: ${FEELPP_ENABLE_PCH}")
  MESSAGE(STATUS "  FEELPP_ENABLE_PCH_APPLICATIONS: ${FEELPP_ENABLE_PCH_APPLICATIONS}")
  MESSAGE(STATUS "            FEELPP_ENABLE_PYTHON: ${FEELPP_ENABLE_PYTHON}")
  MESSAGE(STATUS "   FEELPP_ENABLE_PYTHON_WRAPPING: ${FEELPP_ENABLE_PYTHON_WRAPPING}")
  MESSAGE(STATUS "             FEELPP_ENABLE_CLING: ${FEELPP_ENABLE_CLING}")
  MESSAGE(STATUS " FEELPP_ENABLE_SIMPLE_WEB_SERVER: ${FEELPP_ENABLE_SIMPLE_WEB_SERVER}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "           FEELPP_MESH_MAX_ORDER: ${FEELPP_MESH_MAX_ORDER}")
  MESSAGE(STATUS "       FEELPP_INSTANTIATION_MODE: ${FEELPP_INSTANTIATION_MODE}")
  MESSAGE(STATUS "  FEELPP_INSTANTIATION_ORDER_MAX: ${FEELPP_INSTANTIATION_ORDER_MAX}")
  MESSAGE(STATUS "")
  MESSAGE(STATUS "  FEELPP_CONTRIB_ENABLED_OPTIONS: ${FEELPP_CONTRIB_ENABLED_OPTIONS}")
  MESSAGE(STATUS "          FEELPP_ENABLED_OPTIONS: ${FEELPP_ENABLED_OPTIONS}")
  MESSAGE(STATUS "          FEELPP_ENABLED_MODULES: ${FEELPP_ENABLED_MODULES}")
  MESSAGE(STATUS "         FEELPP_DISABLED_MODULES: ${FEELPP_DISABLED_MODULES}")
  MESSAGE(STATUS "         FEELPP_ENABLED_PROJECTS: ${FEELPP_ENABLED_PROJECTS}")
  MESSAGE(STATUS "                       SCHEDULER: ${FEELPP_SCHEDULER}")
  if (FEELPP_ENABLE_SCHED_OAR )
  message(STATUS "            *** MPIEXEC_PREFLAGS:  ${MPIEXEC_PREFLAGS} ONLY VALID FOR OpenMPI ***")
  endif()

  MESSAGE(STATUS "=================================================================================================")

  string(TOLOWER "${CMAKE_GENERATOR}" cmake_generator_tolower)
  if(cmake_generator_tolower MATCHES "makefile")
    message(STATUS "Some things you can do now with Feel++:")
    MESSAGE(STATUS "===============================================================================================")
    message(STATUS "Command                      |   Description")
    MESSAGE(STATUS "=============================|=================================================================")
    message(STATUS "cd feelpp && make            | Compile the Feel++ library with tools and quickstart applications")
    message(STATUS "cd feelpp && make install    | Compile and install the Feel++ library, tools and apps in ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "cd toolboxes && make         | Compile the Feel++ toolboxes (compiles the feel++ library)")
    message(STATUS "cd toolboxes && make install | Compile and install the Feel++ toolboxes in ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "cd mor && make               | Compile the Feel++ MOR component (compiles the feel++ library)")
    message(STATUS "cd mor && make install       | Compile and install the Feel++ MOR component in ${CMAKE_INSTALL_PREFIX}")  
    message(STATUS "make                         | Compile everything(Feel++ library, toolboxes and mor) ")
    message(STATUS "make install                 | Compile and install everything to ${CMAKE_INSTALL_PREFIX}")
    MESSAGE(STATUS "===============================================================================================")
  endif()

  message(STATUS "")
endfunction(feelpp_print_directive)
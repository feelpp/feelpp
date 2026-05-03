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
  if(FEELPP_ENABLE_VERBOSE_CMAKE)
    message(STATUS "[feelpp] compile definitions for ${varTarget}: ${CD}")
  endif()
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
  MESSAGE(STATUS "Feel++ ${FEELPP_VERSION_STRING} | revision=${FEELPP_REVISION} build=${FEELPP_BUILDID} | c++${FEELPP_STD_CPP}")
  MESSAGE(STATUS "Build | generator=${CMAKE_GENERATOR} type=${CMAKE_BUILD_TYPE} compiler=${CMAKE_CXX_COMPILER}")
  MESSAGE(STATUS "Install | prefix=${CMAKE_INSTALL_PREFIX} data=${FEELPP_DATA_DIR}")
  MESSAGE(STATUS "Features | mpi=${FEELPP_ENABLE_MPI_MODE} python=${FEELPP_ENABLE_PYTHON} wrapping=${FEELPP_ENABLE_PYTHON_WRAPPING} web=${FEELPP_ENABLE_SIMPLE_WEB_SERVER} cling=${FEELPP_ENABLE_CLING}")
  MESSAGE(STATUS "Projects | quickstart=${FEELPP_ENABLE_QUICKSTART} docs=${FEELPP_ENABLE_DOCUMENTATION} doxygen=${FEELPP_ENABLE_DOXYGEN} tests=${FEELPP_ENABLE_TESTS} benchmarks=${FEELPP_ENABLE_BENCHMARKS} research=${FEELPP_ENABLE_RESEARCH}")
  MESSAGE(STATUS "Mesh | max_order=${FEELPP_MESH_MAX_ORDER} instantiation=${FEELPP_INSTANTIATION_MODE} order_max=${FEELPP_INSTANTIATION_ORDER_MAX}")
  MESSAGE(STATUS "Options | contrib=${FEELPP_CONTRIB_ENABLED_OPTIONS}")
  MESSAGE(STATUS "Options | enabled=${FEELPP_ENABLED_OPTIONS}")
  MESSAGE(STATUS "Modules | enabled=${FEELPP_ENABLED_MODULES} disabled=${FEELPP_DISABLED_MODULES}")
  MESSAGE(STATUS "Projects | extra=${FEELPP_ENABLED_PROJECTS} scheduler=${FEELPP_SCHEDULER}")
  if(FEELPP_ENABLE_VERBOSE_CMAKE)
    MESSAGE(STATUS "Verbose | required_flags=${CMAKE_REQUIRED_FLAGS}")
    MESSAGE(STATUS "Verbose | cxx_flags=${CMAKE_CXX_FLAGS}")
    MESSAGE(STATUS "Verbose | cxx_release=${CMAKE_CXX_FLAGS_RELEASE}")
    MESSAGE(STATUS "Verbose | cxx_relwithdebinfo=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    MESSAGE(STATUS "Verbose | cxx_debug=${CMAKE_CXX_FLAGS_DEBUG}")
    MESSAGE(STATUS "Verbose | cxx_debugfull=${CMAKE_CXX_FLAGS_DEBUGFULL}")
    MESSAGE(STATUS "Verbose | cxx_asan=${CMAKE_CXX_FLAGS_ASAN}")
    MESSAGE(STATUS "Verbose | cxx_coverage=${CMAKE_CXX_FLAGS_COVERAGE}")
    MESSAGE(STATUS "Verbose | include_dir=${FEELPP_INCLUDE_DIR}")
    MESSAGE(STATUS "Verbose | machine=${FEELPP_MACHINE_NAME} git=${FEELPP_ENABLE_GIT} nproc=${N} build_nproc=${N2}")
  endif()
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

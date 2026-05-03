###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-28
#
#  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)
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
cmake_minimum_required(VERSION 3.21)
project(testsuite)
set(PROJECT_SHORTNAME "test")

if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "testsuite")
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[testsuite] found Feel++")
  endif()
  SET(FEELPP_DEFAULT_TEST_TIMEOUT "300")
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  feelpp_get_environment()
endif()

enable_testing()
include(CTest)

# Label all tests in this directory tree so they can be run via component targets
set_directory_properties(PROPERTIES LABEL testsuite LABELS "testsuite;testsuite-feelpp")

if ( NOT TARGET testsuite-feelpp )
  add_custom_target(testsuite-feelpp
    #COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L testsuite-feelpp --test-dir ${CMAKE_BINARY_DIR}
    #USES_TERMINAL
  )
endif()

if ( NOT TARGET testsuite )
  add_custom_target(testsuite)
  add_dependencies(testsuite testsuite-feelpp)
endif()

set(FEELPP_HDIV_SMOKE_REGEX
    "feelpp_test_interpolation_RT-np-1\\|feelpp_test_hdiv_oneelt-np-1\\|feelpp_test_hdiv3D_oneelt-np-1\\|feelpp_test_continuity_hdiv-np-1")
set(FEELPP_HDIV_EXTENDED_REGEX
    "feelpp_test_interpolation_BDM-np-1\\|feelpp_test_bdm_static-np-1\\|feelpp_test_rt0_static-np-1\\|feelpp_test_interpolation_Nedelec-np-1\\|feelpp_test_interpolation_Nedelec3d-np-1\\|feelpp_test_hcurl_oneelt-np-1\\|feelpp_test_hcurl3D_oneelt-np-1\\|feelpp_test_continuity_hcurl-np-1\\|feelpp_test_hdiv_runtime_order_parity-np-1\\|feelpp_test_dynamic_order_math_parity-np-1")

set(FEELPP_HDIV_SMOKE_BUDGET_SECONDS 120 CACHE STRING "Expected runtime budget (seconds) for feelpp_hdiv_smoke")
set(FEELPP_HDIV_EXTENDED_BUDGET_SECONDS 900 CACHE STRING "Expected runtime budget (seconds) for feelpp_hdiv_extended")

if ( NOT TARGET feelpp_hdiv_smoke )
  add_custom_target(feelpp_hdiv_smoke
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "${FEELPP_HDIV_SMOKE_REGEX}" --test-dir ${CMAKE_BINARY_DIR}
    USES_TERMINAL
  )
endif()

if ( NOT TARGET feelpp_hdiv_extended )
  add_custom_target(feelpp_hdiv_extended
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R "${FEELPP_HDIV_EXTENDED_REGEX}" --test-dir ${CMAKE_BINARY_DIR}
    USES_TERMINAL
  )
  add_dependencies(feelpp_hdiv_extended feelpp_hdiv_smoke)
endif()

include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/..)

option(FEELPP_DISABLE_COMPILE_FAILURE "disable compilation failure in testsuite" OFF)
if (FEELPP_DISABLE_COMPILE_FAILURE )
  add_definitions(-DFEELPP_DISABLE_COMPILE_FAILURE=1)
endif()
message(STATUS "[feelpp/testsuite] disable compile failure: ${FEELPP_DISABLE_COMPILE_FAILURE}")

find_package(Python3 COMPONENTS Interpreter Development )
message(STATUS "[feelpp/testsuite] Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "[feelpp/testsuite] Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "[feelpp/testsuite] Python3_LIBRARIES: ${Python3_LIBRARIES}")
message(STATUS "[feelpp/testsuite] Python3_VERSION: ${Python3_VERSION}")

foreach( testdir core python math alg mesh poly filters material models interpolation integration discr ts vf leaks opt ls pde fit fmi event viewfactor)
  add_custom_target( test${testdir} )
  add_subdirectory( feel${testdir} )

endforeach()

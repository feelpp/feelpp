
###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-12
#
#  Copyright (C) 2009-2017 Universit� Joseph Fourier (Grenoble I)
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
if(DEFINED PROJECT_NAME)
  set( PYTHON_SUBPROJECT ON )
endif()

project(feelpp-python)
  
set(FEELPP_PYTHON_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(FEELPP_PYTHON_BINARY_DIR ${PROJECT_BINARY_DIR})
SET( CMAKE_MODULE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# with this it is possible to compile applications/MOR from an installed
# version of feel++. This is useful also for staged compilation for docker
list(FIND FEELPP_COMPONENT "mor" COMPONENT_MOR)
list(FIND FEELPP_COMPONENT "python" COMPONENT_PYTHON )

if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR COMPONENT_PYTHON GREATER -1 )
  if ( COMPONENT_MOR GREATER -1 )
    find_package(Feel++ COMPONENTS Toolboxes PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
    set(FEELPP_MOR_FOUND TRUE)
  else()
    find_package(Feel++ COMPONENTS Toolboxes MOR PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules)
  endif()
    if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "[python] Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[python] found Feel++ library")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(WARNING "[python] Feel++ was not installed with toolboxes component")
  else()
    message(STATUS "[python] found Feel++ toolboxes")
  endif()
    if (NOT FEELPP_MOR_FOUND)
    message(WARNING "[python] Feel++ was not installed with mor component")
  else()
    message(STATUS "[python] found Feel++ mor")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  feelpp_get_environment()
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
  set(FEELPP_FOUND ON)
  set(FEELPP_TOOLBOXES_FOUND ${FEELPP_ENABLE_TOOLBOXES})
  set(FEELPP_MOR_FOUND ${FEELPP_ENABLE_MOR})
endif()

find_package(Python3 COMPONENTS Interpreter Development)
message(STATUS "[feelpp/python] Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "[feelpp/python] Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "[feelpp/python] Python3_LIBRARIES: ${Python3_LIBRARIES}")
message(STATUS "[feelpp/python] Python3_VERSION: ${Python3_VERSION}")

if ( FEELPP_HAS_PYTHON )
  message(STATUS "[feelpp-python] python module path:${FEELPP_PYTHON_MODULE_PATH}")
  find_package(MPI4PY REQUIRED)
  if ( MPI4PY_FOUND )
    message(STATUS "[feelpp-python] mpi4py header ${MPI4PY_INCLUDE_DIR}")
  endif()
  find_package(PETSC4PY)
  if ( PETSC4PY_FOUND )
    message(STATUS "[feelpp-python] petsc4py header ${PETSC4PY_INCLUDE_DIR}")
  endif()
  find_package(SLEPC4PY)
  if ( SLEPC4PY_FOUND )
    message(STATUS "[feelpp-python] slepc4py header ${SLEPC4PY_INCLUDE_DIR}")
  endif()
endif()

#
# add a pybind11 feelpp module
# FEELPP_PYTHON_MODULE_PATH must be defined !
#
macro(feelpp_add_pymodule)
 PARSE_ARGUMENTS(FEELPP_PYMODULE
    "NAME;SRCS;DESTINATION;LINK_LIBRARIES;INCLUDES"
    ""
    ${ARGN}
    )
  CAR(FEELPP_PYMODULE_NAME ${FEELPP_PYMODULE_DEFAULT_ARGS})
  message(STATUS "[pyfeelpp] add pymodule ${FEELPP_PYMODULE_NAME}")
  pybind11_add_module(_${FEELPP_PYMODULE_NAME}  ${FEELPP_PYMODULE_SRCS}  )
  target_include_directories(_${FEELPP_PYMODULE_NAME} PRIVATE ${PYTHON_INCLUDE_DIRS} ${MPI4PY_INCLUDE_DIR} ${PETSC4PY_INCLUDE_DIR} ${SLEPC4PY_INCLUDE_DIR})
  target_link_libraries( _${FEELPP_PYMODULE_NAME} PUBLIC Feelpp::feelpp ${FEELPP_PYMODULE_LINK_LIBRARIES} )
  install(TARGETS _${FEELPP_PYMODULE_NAME} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/${FEELPP_PYMODULE_DESTINATION})
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py )
    add_custom_command(
           TARGET _${FEELPP_PYMODULE_NAME} POST_BUILD
           COMMAND ${CMAKE_COMMAND} -E copy
                   ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
                   ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
  endif()
endmacro(feelpp_add_pymodule)

option(FEELPP_ENABLE_PYTHON3_VENV ON)
if (FEELPP_ENABLE_PYTHON3_VENV)
  set(REQUIREMENTS ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt)
  set(VENV ${CMAKE_CURRENT_BINARY_DIR}/venv)
  set(BIN_DIR ${VENV}/bin)
  set(PYTHON ${BIN_DIR}/python)
  
  add_custom_command(
    OUTPUT ${PYTHON}
    DEPENDS venv requirements.txt
    COMMAND ${Python3_EXECUTABLE} -m venv ${VENV}
    COMMAND ${BIN_DIR}/pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt 
    DEPENDS ${REQUIREMENTS}
  )
  add_custom_target(venv DEPENDS ${PYTHON})
  #add_dependencies(test venv)
endif(FEELPP_ENABLE_PYTHON3_VENV)

#include(feelpp.module.pyfeelpp)
if ( FEELPP_FOUND )
add_subdirectory(pyfeelpp)
add_subdirectory(notebooks)
endif()

# TOOLBOXES needs to be defined for pyfeelpp-toolboxes 
if (FEELPP_TOOLBOXES_FOUND )
add_subdirectory( pyfeelpp-toolboxes )
endif()

if ( FEELPP_MOR_FOUND ) 
add_subdirectory(pyfeelpp-mor)
endif()
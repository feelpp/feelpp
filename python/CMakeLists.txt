
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
message(WARNING "[python] legacy monolithic Python entry point kept for compatibility; prefer feelpp/python, toolboxes/python, and mor/python")
  
set(FEELPP_PYTHON_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(FEELPP_PYTHON_BINARY_DIR ${PROJECT_BINARY_DIR})
SET( CMAKE_MODULE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
get_filename_component(FEELPP_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

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

find_package(Python3 COMPONENTS Interpreter Development Development.Module)
message(STATUS "[feelpp/python] Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "[feelpp/python] Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "[feelpp/python] Python3_LIBRARIES: ${Python3_LIBRARIES}")
message(STATUS "[feelpp/python] Python3_VERSION: ${Python3_VERSION}")
include(FeelppPythonNormalize)

# If the project is configured to use an external pybind11, make sure the
# expected imported targets exist before any pybind11_add_module() calls.
if(DEFINED FEELPP_USE_EXTERNAL_PYBIND11 AND FEELPP_USE_EXTERNAL_PYBIND11)
  find_package(pybind11 CONFIG REQUIRED)
  if(NOT TARGET pybind11::headers)
    if(TARGET pybind11::pybind11_headers)
      add_library(pybind11::headers ALIAS pybind11::pybind11_headers)
    elseif(TARGET pybind11::module)
      add_library(pybind11::headers ALIAS pybind11::module)
    elseif(TARGET pybind11::pybind11)
      add_library(pybind11::headers ALIAS pybind11::pybind11)
    elseif(TARGET pybind11)
      add_library(pybind11::headers ALIAS pybind11)
    endif()
  endif()
endif()


feelpp_python_normalize_module_path()


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

include("${FEELPP_ROOT_SOURCE_DIR}/feelpp/cmake/modules/FeelppPythonSupport.cmake")

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
add_subdirectory("${FEELPP_ROOT_SOURCE_DIR}/feelpp/python" "${CMAKE_BINARY_DIR}/feelpp/python")
add_subdirectory(notebooks)
endif()

# TOOLBOXES needs to be defined for pyfeelpp-toolboxes 
if (FEELPP_TOOLBOXES_FOUND )
add_subdirectory("${FEELPP_ROOT_SOURCE_DIR}/toolboxes/python" "${CMAKE_BINARY_DIR}/toolboxes/python")
endif()

if ( FEELPP_MOR_FOUND ) 
add_subdirectory("${FEELPP_ROOT_SOURCE_DIR}/mor/python" "${CMAKE_BINARY_DIR}/mor/python")
endif()

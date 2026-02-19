
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


# When building the Python component standalone against an existing Feel++
# installation, FEELPP_PYTHON_MODULE_PATH may come pre-baked for a different
# prefix (e.g. an absolute "/usr/local/..." path or a relative "local/..."
# path). For packaging, honor an explicit PYTHON_SITE_PACKAGES override and
# otherwise recompute the module path for the current CMAKE_INSTALL_PREFIX
# using Python's sysconfig.
if (DEFINED PYTHON_SITE_PACKAGES)
  # Respect an explicit override, but guard against absolute paths under a different
  # prefix (e.g. /usr/local/...) which can come from an installed Feel++ config.
  set(_REQ_PY_SITE "${PYTHON_SITE_PACKAGES}")
  set(_USE_REQ_PATH ON)
  if (IS_ABSOLUTE "${_REQ_PY_SITE}")
    string(FIND "${_REQ_PY_SITE}" "${CMAKE_INSTALL_PREFIX}/" _req_prefix_pos)
    if (_req_prefix_pos EQUAL -1)
      set(_USE_REQ_PATH OFF)
    else()
      set(_LOCAL_PREFIX "${CMAKE_INSTALL_PREFIX}/local/")
      string(FIND "${_REQ_PY_SITE}" "${_LOCAL_PREFIX}" _req_local_pos)
      if (_req_local_pos EQUAL 0)
        set(_USE_REQ_PATH OFF)
      endif()
    endif()
  endif()
  if (_USE_REQ_PATH)
    # Ensure both normal and cache variables are updated to avoid scope shadowing
    set(FEELPP_PYTHON_MODULE_PATH "${_REQ_PY_SITE}")
    set(FEELPP_PYTHON_MODULE_PATH "${_REQ_PY_SITE}" CACHE STRING "Python module install path" FORCE)
  else()
    # Recompute from sysconfig for the current prefix
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "\nimport sys, sysconfig, site, os\nbase='${CMAKE_INSTALL_PREFIX}'\ncands=[]\ntry:\n    cands+=site.getsitepackages()\nexcept Exception:\n    pass\ncands.append(sysconfig.get_path('platlib', vars={'base': base, 'platbase': base}))\n# Prefer Debian-style '/usr/lib/python3/dist-packages' under current base\npreferred=[base + '/lib/python3/dist-packages'] \nchoice=None\nfor p in cands:\n    if not p: continue\n    if any(p.startswith(pr) for pr in preferred):\n        choice=p; break\n# Next, prefer any path under the current base and not under base+'/local'\nif choice is None:\n    for p in cands:\n        if p.startswith(base+os.sep) and not p.startswith(base+os.sep+'local'+os.sep):\n            choice=p; break\n# Fallback to first candidate\nif choice is None:\n    choice=cands[0]\nprint(choice)\n"
      OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT "${_ABS_PYTHON_MODULE_PATH}" STREQUAL "")
      get_filename_component(_ABS_PYTHON_MODULE_PATH "${_ABS_PYTHON_MODULE_PATH}" ABSOLUTE)
      file(RELATIVE_PATH FEELPP_PYTHON_MODULE_PATH "${CMAKE_INSTALL_PREFIX}" "${_ABS_PYTHON_MODULE_PATH}")
      set(FEELPP_PYTHON_MODULE_PATH "${FEELPP_PYTHON_MODULE_PATH}" CACHE STRING "Python module install path" FORCE)
    endif()
  endif()
elseif (DEFINED FEELPP_PYTHON_MODULE_PATH)
  set(_RECOMPUTE_PATH OFF)
  # Case 1: relative path starting with "local/"
  if (FEELPP_PYTHON_MODULE_PATH MATCHES "^local/")
    set(_RECOMPUTE_PATH ON)
  endif()
  # Case 2: absolute path under a different prefix than CMAKE_INSTALL_PREFIX
  if (IS_ABSOLUTE "${FEELPP_PYTHON_MODULE_PATH}")
    # Recompute if it does not live under the current prefix OR
    # if it points to the "local" subtree (e.g. /usr/local/...)
    string(FIND "${FEELPP_PYTHON_MODULE_PATH}" "${CMAKE_INSTALL_PREFIX}/" _prefix_pos)
    if (_prefix_pos EQUAL -1)
      set(_RECOMPUTE_PATH ON)
    else()
      # Also recompute when prefix is "/usr" but path is under "/usr/local"
      set(_LOCAL_PREFIX "${CMAKE_INSTALL_PREFIX}/local/")
      string(FIND "${FEELPP_PYTHON_MODULE_PATH}" "${_LOCAL_PREFIX}" _local_pos)
      if (_local_pos EQUAL 0)
        set(_RECOMPUTE_PATH ON)
      endif()
    endif()
  endif()
  if (_RECOMPUTE_PATH)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "\nimport sys, sysconfig, site, os\nbase='${CMAKE_INSTALL_PREFIX}'\ncands=[]\ntry:\n    cands+=site.getsitepackages()\nexcept Exception:\n    pass\ncands.append(sysconfig.get_path('platlib', vars={'base': base, 'platbase': base}))\npreferred=[base + '/lib/python3/dist-packages'] \nchoice=None\nfor p in cands:\n    if not p: continue\n    if any(p.startswith(pr) for pr in preferred):\n        choice=p; break\nif choice is None:\n    for p in cands:\n        if p.startswith(base+os.sep) and not p.startswith(base+os.sep+'local'+os.sep):\n            choice=p; break\nif choice is None:\n    choice=cands[0]\nprint(choice)\n"
      OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT "${_ABS_PYTHON_MODULE_PATH}" STREQUAL "")
      get_filename_component(_ABS_PYTHON_MODULE_PATH "${_ABS_PYTHON_MODULE_PATH}" ABSOLUTE)
      file(RELATIVE_PATH FEELPP_PYTHON_MODULE_PATH "${CMAKE_INSTALL_PREFIX}" "${_ABS_PYTHON_MODULE_PATH}")
      set(FEELPP_PYTHON_MODULE_PATH "${FEELPP_PYTHON_MODULE_PATH}" CACHE STRING "Python module install path" FORCE)
    endif()
  endif()
else()
  # Neither override nor pre-baked value given: compute from sysconfig for current prefix
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "\nimport sys, sysconfig, site, os\nbase='${CMAKE_INSTALL_PREFIX}'\ncands=[]\ntry:\n    cands+=site.getsitepackages()\nexcept Exception:\n    pass\ncands.append(sysconfig.get_path('platlib', vars={'base': base, 'platbase': base}))\npreferred=[base + '/lib/python3/dist-packages'] \nchoice=None\nfor p in cands:\n    if not p: continue\n    if any(p.startswith(pr) for pr in preferred):\n        choice=p; break\nif choice is None:\n    for p in cands:\n        if p.startswith(base+os.sep) and not p.startswith(base+os.sep+'local'+os.sep):\n            choice=p; break\nif choice is None:\n    choice=cands[0]\nprint(choice)\n"
    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT "${_ABS_PYTHON_MODULE_PATH}" STREQUAL "")
    get_filename_component(_ABS_PYTHON_MODULE_PATH "${_ABS_PYTHON_MODULE_PATH}" ABSOLUTE)
    file(RELATIVE_PATH FEELPP_PYTHON_MODULE_PATH "${CMAKE_INSTALL_PREFIX}" "${_ABS_PYTHON_MODULE_PATH}")
    set(FEELPP_PYTHON_MODULE_PATH "${FEELPP_PYTHON_MODULE_PATH}" CACHE STRING "Python module install path" FORCE)
  endif()
endif()


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

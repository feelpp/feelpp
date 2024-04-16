###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 19 Jun 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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
cmake_minimum_required(VERSION 3.13)
project(pyfeelpp-toolboxes)
if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "toolboxes/pyfeelpp-toolboxes" )
  find_package(Feel++ COMPONENTS Toolboxes PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "[pyfeelpp-toolboxes] Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[pyfeelpp-toolboxes] found Feel++ library")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(WARNING "[pyfeelpp-toolboxes] Feel++ was not installed with toolboxes component")
  else()
    message(STATUS "[pyfeelpp-toolboxes] found Feel++ toolboxes")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
endif()
#
# setup pyfeelpp-toolboxes
#
OPTION(FEELPP_ENABLE_PYFEELPP_TOOLBOXES "Enable PyFeel++ for toolboxes" ON  )
if (FEELPP_ENABLE_PYFEELPP_TOOLBOXES)
  message(STATUS "Enabling PyFeel++ Toolboxes" )
  add_subdirectory(feelpp)
  add_subdirectory(tests)
  # additional cmake args can be necessary for paraview
  # set(PYFEELPP_SETUP_HAS_PARAVIEW_CMAKE_ARGS 0)
  # set(PYFEELPP_SETUP_PARAVIEW_CMAKE_ARGS "[]")
  # if ( FEELPP_HAS_PARAVIEW AND PARAVIEW_ENABLE_QT_SUPPORT )
  #   if (PARAVIEW_QT_VERSION VERSION_GREATER "4")
  #     set(PYFEELPP_SETUP_HAS_PARAVIEW_CMAKE_ARGS 1)
  #     set(PYFEELPP_SETUP_PARAVIEW_CMAKE_ARGS "['-DQt${PARAVIEW_QT_VERSION}_DIR=${Qt${PARAVIEW_QT_VERSION}_DIR}']")
  #   else()
  #     # TODO QT4
  #   endif()
  # endif()
# 
  # set(FEELPP_TOOLBOXES_DIR ${CMAKE_INSTALL_PREFIX})
# 
  # set(SETUP_TOOLBOXES_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup-toolboxes.py")
  # set(SETUP_TOOLBOXES_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup-toolboxes.py")
  # set(PYFEELPP_TOOLBOXES_PYTHONPATH $ENV{PYTHONPATH})
  # configure_file(${SETUP_TOOLBOXES_PY_IN} ${SETUP_TOOLBOXES_PY})
  # set(PYFEELPP_TOOLBOXES_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/)
  # install(DIRECTORY DESTINATION ${PYFEELPP_TOOLBOXES_INSTALL_DIR})
  # install(
  #   CODE "set(ENV{PYTHONPATH} \"${PYFEELPP_TOOLBOXES_INSTALL_DIR}:$ENV{PYTHONPATH}\")"
  #   CODE "
  #     execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_TOOLBOXES_PY} install --prefix ${CMAKE_INSTALL_PREFIX} RESULT_VARIABLE ret )
  #     if(ret EQUAL \"1\")
  #       message( FATAL_ERROR \"pyfeelpp-toolboxes setup.py bad exit status\")
  #     endif()
  #   "
  # )
endif()

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
cmake_minimum_required(VERSION 3.3)
project(pyfeelpp)

if ( ${CMAKE_SOURCE_DIR}/ STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
else()
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
endif()

# additional cmake args can be necessary for paraview
set(PYFEELPP_SETUP_HAS_PARAVIEW_CMAKE_ARGS 0)
set(PYFEELPP_SETUP_PARAVIEW_CMAKE_ARGS "[]")
if ( FEELPP_HAS_PARAVIEW AND PARAVIEW_ENABLE_QT_SUPPORT )
  if (PARAVIEW_QT_VERSION VERSION_GREATER "4")
    set(PYFEELPP_SETUP_HAS_PARAVIEW_CMAKE_ARGS 1)
    set(PYFEELPP_SETUP_PARAVIEW_CMAKE_ARGS "['-DQt${PARAVIEW_QT_VERSION}_DIR=${Qt${PARAVIEW_QT_VERSION}_DIR}']")
  else()
    # TODO QT4
  endif()
endif()

enable_testing()

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  list(APPEND CMAKE_MODULE_PATH "$ENV{FEELPP_DIR}/share/cmake/pybind11")
  list(APPEND CMAKE_MODULE_PATH "/usr/local/share/cmake/pybind11")
  include(pybind11Tools)
endif()

#
# setup for pyfeelpp-libfeelpp
#
# OPTION(FEELPP_ENABLE_PYFEELPP_LIBFEELPP "Enable PyFeel++ for libfeelpp" ON  )
# if (FEELPP_ENABLE_PYFEELPP_LIBFEELPP)
#   set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
#   set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
#   configure_file(${SETUP_PY_IN} ${SETUP_PY})
#   set(PYFEELPP_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/)
#   install(DIRECTORY DESTINATION ${PYFEELPP_INSTALL_DIR})
#   install(
#     CODE "set(ENV{PYTHONPATH} \"${PYFEELPP_INSTALL_DIR}:$ENV{PYTHONPATH}\")"
#     CODE "
#     execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix ${CMAKE_INSTALL_PREFIX} RESULT_VARIABLE ret )
#     if(ret EQUAL \"1\")
#        message( FATAL_ERROR \"pyfeelpp setup.py bad exit status\")
#     endif()
#     "
#     )
# endif()

#  if (WIN32)
#    set (PYTHON_SITE_PACKAGES Lib/site-packages/feelpp)
#  else ()
#    set (PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/feelpp)
#  endif ()

add_subdirectory(feelpp)
add_subdirectory(tests)
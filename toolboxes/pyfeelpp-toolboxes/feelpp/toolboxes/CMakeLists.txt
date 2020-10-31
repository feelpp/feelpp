###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 06 Aug 2018
#
#  Copyright (C) 2018 Feel++ Consortium
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

add_custom_target(pyfeelpptoolboxes)

if (DEFINED PYTHON_SITE_PACKAGES)
  set (FEELPP_PYTHON_MODULE_PATH ${PYTHON_SITE_PACKAGES})
else ()
  execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
                    OUTPUT_VARIABLE _ABS_PYTHON_MODULE_PATH
                    RESULT_VARIABLE _PYTHON_pythonlib_result
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if (_PYTHON_pythonlib_result)
    message (SEND_ERROR "Could not run ${PYTHON_EXECUTABLE}")
  endif ()

  get_filename_component (_ABS_PYTHON_MODULE_PATH ${_ABS_PYTHON_MODULE_PATH} ABSOLUTE)
  file (RELATIVE_PATH FEELPP_PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_MODULE_PATH})
endif ()
set (FEELPP_PYTHON${PYTHON_VERSION_MAJOR}_MODULE_PATH ${FEELPP_PYTHON_MODULE_PATH})
message(STATUS "[feelpp-toolboxes] python module path: ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes")

#add_dependencies(pyfeelpptoolboxes pyfeelpp)

pybind11_add_module(_modelcore  modelcore.cpp)
target_link_libraries( _modelcore PRIVATE Feelpp::feelpp_modelcore )

pybind11_add_module(_modelmesh  modelmesh.cpp)
target_link_libraries( _modelmesh PRIVATE Feelpp::feelpp_modelmesh Feelpp::feelpp_modelmeshale)


install(TARGETS _modelcore _modelmesh DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes/)
install(FILES  __init__.py DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/toolboxes/)



add_dependencies(pyfeelpptoolboxes _modelcore _modelmesh)

if ( 0 )
if ( FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( fluid )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solid )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_FSI )
  add_subdirectory( fsi )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_HDG AND ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )
  add_subdirectory( hdg )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_ADVECTION )
  add_subdirectory( advection )
endif()
if( FEELPP_TOOLBOXES_ENABLE_LEVELSET )
  add_subdirectory( levelset )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MULTIFLUID )
  add_subdirectory( multifluid )
endif()
  # Heat Toolbox

if ( FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( heatfluid )
endif()


if ( FEELPP_TOOLBOXES_ENABLE_MESHALE )
  add_subdirectory( alemesh )
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MAXWELL )
  add_subdirectory( maxwell )
endif()
else()
if ( FEELPP_TOOLBOXES_ENABLE_HEAT )
  add_subdirectory( heat )
  add_dependencies(pyfeelpptoolboxes _heat)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_ELECTRIC )
  add_subdirectory( electric )
  add_dependencies(pyfeelpptoolboxes _electric)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC )
  add_subdirectory( thermoelectric )
  add_dependencies(pyfeelpptoolboxes _thermoelectric)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solid )
  add_dependencies(pyfeelpptoolboxes _solid)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( fluid )
  add_dependencies(pyfeelpptoolboxes _fluid)
endif()
endif()

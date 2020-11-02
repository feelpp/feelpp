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

feelpp_add_pymodule(modelcore SRCS modelcore.cpp DESTINATION feelpp/toolboxes/ LINK_LIBRARIES Feelpp::feelpp_modelcore )
feelpp_add_pymodule(modelmesh SRCS modelmesh.cpp DESTINATION feelpp/toolboxes/ LINK_LIBRARIES Feelpp::feelpp_modelmesh Feelpp::feelpp_modelmeshale)

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

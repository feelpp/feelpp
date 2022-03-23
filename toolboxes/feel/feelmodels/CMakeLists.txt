###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 23 Feb 2015
#
#  Copyright (C) 2015-present Feel++ Consortium
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

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp.models.macros.cmake)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp.models.macros.cmake DESTINATION share/feelpp/toolboxes/cmake COMPONENT Devel)

add_subdirectory(cmake)
add_subdirectory(modelcore)
add_subdirectory(modelmesh)
add_subdirectory(modelvf)
add_subdirectory(modelmaterials)


add_custom_target(feelpp-models-mesh)
add_dependencies(feelpp-models-mesh feelpp_modelmesh)

add_custom_target(feelpp-models)
#add_dependencies(feelpp-models feelpp_model_heat feelpp_model_solidmechanics feelpp_model_fluid feelpp_model_fsi)
add_dependencies(feelpp-models feelpp-models-mesh)

add_custom_target(feelpp-models-common)
add_dependencies(feelpp-models-common feelpp_modelcore feelpp_modelmesh)
if( FEELPP_MODELS_HAS_MESHALE )
  add_dependencies(feelpp-models-common feelpp_modelmeshale)
endif()

# generate models lib
if ( FEELPP_TOOLBOXES_ENABLE_HEAT )
  add_subdirectory(heat)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_ELECTRIC )
  add_subdirectory(electric)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS )
  add_subdirectory(solid)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_CFPDE )
  add_subdirectory(coefficientformpdes)
else()
  set(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS OFF)
  set(FEELPP_TOOLBOXES_ENABLE_FSI OFF)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( body )
  add_subdirectory(fluid)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_FSI )
  add_subdirectory(fsi)
endif()

if ( FEELPP_TOOLBOXES_ENABLE_ADVECTION )
  add_subdirectory(advection)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_LEVELSET )
  add_subdirectory(levelset)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MULTIFLUID )
  add_subdirectory(multifluid)
endif() 
if ( FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC )
  add_subdirectory(thermoelectric)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory(heatfluid)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_HDG  AND ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )
  add_subdirectory(hdg)
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MAXWELL )
  add_subdirectory( maxwell )
endif()


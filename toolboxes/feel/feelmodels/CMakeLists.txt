###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 23 Feb 2015
#
#  Copyright (C) 2015 Feel++ Consortium
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

SET(FEELPP_MODELS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "") # TODO to remove
SET(FEELPP_MODELS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "") # TODO to remove

add_subdirectory(modelcore)
add_subdirectory(modelalg)
add_subdirectory(modelmesh)
add_subdirectory(modelvf)


add_custom_target(feelpp-models-algebra)
add_dependencies(feelpp-models-algebra feelpp_modelalg)

add_custom_target(feelpp-models-mesh)
add_dependencies(feelpp-models-mesh feelpp_modelmesh)

add_custom_target(feelpp-models)
#add_dependencies(feelpp-models feelpp_model_thermodynamics feelpp_model_solidmechanics feelpp_model_fluid feelpp_model_fsi)
add_dependencies(feelpp-models feelpp-models-algebra feelpp-models-mesh)

add_custom_target(feelpp-models-common)
add_dependencies(feelpp-models-common feelpp_modelcore feelpp_modelalg feelpp_modelmesh)
if( FEELPP_MODELS_HAS_MESHALE )
  add_dependencies(feelpp-models-common feelpp_modelmeshale)
endif()

# generate models lib
if ( FEELPP_MODELS_ENABLE_THERMODYNAMICS )
  add_subdirectory(thermodyn)
endif()
if ( FEELPP_MODELS_ENABLE_SOLIDMECHANICS )
  add_subdirectory(solid)
endif()
if ( FEELPP_MODELS_ENABLE_FLUIDMECHANICS )
  add_subdirectory(fluid)
endif()
if ( FEELPP_MODELS_ENABLE_FSI )
  add_subdirectory(fsi)
endif()
if ( 1 ) # enable
if ( FEELPP_MODELS_ENABLE_ADVECTION )
  add_subdirectory(advection)
endif()
if ( FEELPP_MODELS_ENABLE_LEVELSET )
  add_subdirectory(levelset)
endif()
if ( FEELPP_MODELS_ENABLE_MULTIFLUID )
  add_subdirectory(multifluid)
endif()
endif() 
if ( FEELPP_MODELS_ENABLE_THERMOELECTRIC )
  add_subdirectory(thermoelectric)
endif()
if ( FEELPP_MODELS_ENABLE_HDG )
  add_subdirectory(hdg)
endif()

set( FEELPP_TOOLBOXES_CONFIG_FILE  ${FEELPP_TOOLBOXES_BINARY_DIR}/feel/feelmodels/cmake/Feel++-ToolboxesConfig.cmake)
configure_file(${FEELPP_TOOLBOXES_SOURCE_DIR}/feel/feelmodels/cmake/Feel++-ToolboxesConfig.cmake.in ${FEELPP_TOOLBOXES_CONFIG_FILE} @ONLY)
install(FILES ${FEELPP_TOOLBOXES_CONFIG_FILE} DESTINATION share/feelpp/toolboxes/cmake COMPONENT Devel)
install(EXPORT feelpp-toolboxes-targets DESTINATION share/feelpp/toolboxes/cmake )

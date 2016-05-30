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
SET( SRC
  modelproperties.cpp modelparameters.cpp modelmaterials.cpp modelpostprocess.cpp modelfunctions.cpp
  )

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp solid/*.hpp fluid/*.hpp )
feel_append_src(feelmodels "${SRC};${HDR}")

add_subdirectory(solid_oldversion)

SET(FEELPP_MODELS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_MODELS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")


option(FEELPP_MODELS_ENABLE_MESHALE "Enable/Disable mesh ALE" ON)
if ( FEELPP_MODELS_ENABLE_MESHALE )
  set( FEELPP_MODELS_HAS_MESHALE 1)
endif()
set(FEELPP_MODELS_HAS_MESHALE_HARMONICEXTENSION 1)
#set(FEELPP_MODELS_HAS_MESHALE_WINSLOW 1)

option(FEELPP_MODELS_ENABLE_THERMODYNAMICS "Enable/Disable ThermoDynamics" ON)
option(FEELPP_MODELS_ENABLE_FLUIDMECHANICS "Enable/Disable FluidMechanics" ON)
option(FEELPP_MODELS_ENABLE_SOLIDMECHANICS "Enable/Disable SolidMechanics" ON)
option(FEELPP_MODELS_ENABLE_FSI "Enable/Disable FSI" ON)
option(FEELPP_MODELS_ENABLE_ADVECTION "Enable/Disable Advection" ON)
option(FEELPP_MODELS_ENABLE_LEVELSET "Enable/Disable LevelSet" ON)


add_subdirectory(modelcore)
add_subdirectory(modelalg)
add_subdirectory(modelmesh)


add_custom_target(feelpp-models-algebra)
add_dependencies(feelpp-models-algebra feelpp_modelalg)

add_custom_target(feelpp-models-mesh)
add_dependencies(feelpp-models-mesh feelpp_modelmesh)

add_custom_target(feelpp-models)
#add_dependencies(feelpp-models feelpp_model_thermodynamics feelpp_model_solidmechanics feelpp_model_fluid feelpp_model_fsi)
add_dependencies(feelpp-models feelpp-models-algebra feelpp-models-mesh)


set(FEELPP_MODELS_INSTALL_LIBBASE_DEPENDS  feelpp_modelcore feelpp_modelmesh )
if( FEELPP_MODELS_HAS_MESHALE )
  set(FEELPP_MODELS_INSTALL_LIBBASE_DEPENDS feelpp_modelmeshale ${FEELPP_MODELS_INSTALL_LIBBASE_DEPENDS} )
endif()
add_custom_target(install-feelpp-models-common
  DEPENDS ${FEELPP_MODELS_INSTALL_LIBBASE_DEPENDS}
  COMMAND
      "${CMAKE_COMMAND}" -DCMAKE_INSTALL_COMPONENT=LibsFeelppModelsCommon
      -P "${CMAKE_BINARY_DIR}/cmake_install.cmake"
)


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
if ( FEELPP_MODELS_ENABLE_ADVECTION )
  add_subdirectory(advection)
endif()
if ( FEELPP_MODELS_ENABLE_LEVELSET )
  add_subdirectory(levelset)
endif()

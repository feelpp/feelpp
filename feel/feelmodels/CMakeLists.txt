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
  modelproperties.cpp modelparameters.cpp modelmaterials.cpp modelpostprocess.cpp

  #solid/solidmechanics.cpp
  )

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp solid/*.hpp fluid/*.hpp )
feel_append_src(feelmodels "${SRC};${HDR}")

add_subdirectory(solid_oldversion)

SET(FEELPP_MODELS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(FEELPP_MODELS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")

option(FEELMODELS_ENABLE_MESHALE "Enable/Disable mesh ALE" ON)
if ( FEELMODELS_ENABLE_MESHALE )
  set( FEELPP_MODELS_HAS_MESHALE 1)
endif()
set(FEELPP_MODELS_HAS_MESHALE_HARMONICEXTENSION 1)
#set(FEELPP_MODELS_HAS_MESHALE_WINSLOW 1)


add_subdirectory(modelcore)
add_subdirectory(modelalg)
add_subdirectory(modelmesh)


# generate models lib
add_subdirectory(thermodyn)
add_subdirectory(solid)
add_subdirectory(fluid)
add_subdirectory(fsi)

add_custom_target(feelpp_models)
#add_dependencies(feelpp_models feelpp_model_thermodyn feelpp_model_solid feelpp_model_fluid feelpp_model_fsi)
add_dependencies(feelpp_models feelpp_model_fluid)

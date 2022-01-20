###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 05 Jan 2017
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

project(toolboxes)
set(PROJECT_SHORTNAME "toolbox")

set(FEELPP_TOOLBOXES_SHARED_VERSION 1.0)
set(FEELPP_TOOLBOXES_SHARED_SOVERSION 1.0)

# with this it is possible to compile applications/models from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR FEELPP_COMPONENT STREQUAL "toolboxes")
  find_package(Feel++ PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  else()
    message(STATUS "[feelpp-toolboxes] found Feel++ library")
  endif()
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "install prefix" FORCE)
    # SKIP RPATH if we set CMAKE_INSTALL_PREFIX
    SET(CMAKE_SKIP_RPATH ON)
  endif()
else()
  feelpp_get_environment()
  # usefull for pyfeelpp-toolboxes
  set(FEELPP_DIR ${CMAKE_INSTALL_PREFIX})
endif()

if ( FEELPP_HAS_PYTHON )
  message(STATUS "[feelpp-toolboxes] python module path:${FEELPP_PYTHON_MODULE_PATH}")
endif()


option( FEELPP_ENABLE_ASCIIDOCTOR "Enable AsciiDoctor Support" ${FEELPP_ENABLE_PACKAGE_DEFAULT_OPTION} )

if ( FEELPP_ENABLE_ASCIIDOCTOR )
  include( feelpp.adoc )
endif()

option(FEELPP_TOOLBOXES_ENABLE_MESHALE "Enable/Disable mesh ALE" ON)
if ( FEELPP_TOOLBOXES_ENABLE_MESHALE )
  set(FEELPP_MODELS_HAS_MESHALE 1)
  set(FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION 1)
  set(FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW 1)
endif()

option(FEELPP_TOOLBOXES_ENABLE "Enable/Disable Toolboxes" ON)
option(FEELPP_TOOLBOXES_ENABLE_HEAT "Enable/Disable Heat" ON)
option(FEELPP_TOOLBOXES_ENABLE_CFPDE "Enable/Disable CFPDE" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_ELECTRIC "Enable/Disable Electric" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC "Enable/Disable ThermoElectric" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS "Enable/Disable FluidMechanics" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS "Enable/Disable SolidMechanics" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_FSI "Enable/Disable FSI" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_ADVECTION "Enable/Disable Advection" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_LEVELSET "Enable/Disable LevelSet" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_MULTIFLUID "Enable/Disable MultiFluid" OFF ) # ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_HDG "Enable/Disable HDG" ${FEELPP_TOOLBOXES_ENABLE})
option(FEELPP_TOOLBOXES_ENABLE_MAXWELL "Enable/Disable Magnetostatic" ${FEELPP_TOOLBOXES_ENABLE})
#option(FEELPP_ENABLE_PCH_MODELS "Enable precompiled headers for models (pch)" OFF)

set(FEELPP_TOOLBOXES_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(FEELPP_TOOLBOXES_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

#include_directories(${CMAKE_SOURCE_DIR}/feelpp)
#include_directories(${FEELPP_TOOLBOXES_SOURCE_DIR})
#include_directories(${FEELPP_TOOLBOXES_BINARY_DIR})

INCLUDE(CTest)
ENABLE_TESTING()

add_subdirectory(cmake)
add_subdirectory( feel/feelmodels )

list(APPEND TOOLBOXES_CORE "modelalg" "modelmesh" "modelcore")

if ( FEELPP_TOOLBOXES_ENABLE_CFPDE )
  add_subdirectory( coefficientformpdes )
  list(APPEND TOOLBOXES "coefficientformpdes")
endif()

if ( FEELPP_TOOLBOXES_ENABLE_CFPDE AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( fluid )
  list(APPEND TOOLBOXES "fluid")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solid )
  list(APPEND TOOLBOXES "solid")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_CFPDE AND FEELPP_TOOLBOXES_ENABLE_FSI )
  add_subdirectory( fsi )
  list(APPEND TOOLBOXES "fsi")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_HDG AND ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )
  add_subdirectory( hdg )
  list(APPEND TOOLBOXES "hdg")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_ADVECTION )
  add_subdirectory( advection )
  list(APPEND TOOLBOXES "advection")
endif()
if( FEELPP_TOOLBOXES_ENABLE_LEVELSET )
  add_subdirectory( levelset )
  list(APPEND TOOLBOXES "levelset")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MULTIFLUID )
  add_subdirectory( multifluid )
  list(APPEND TOOLBOXES "multifluid")
endif()
  # Heat Toolbox
if ( FEELPP_TOOLBOXES_ENABLE_HEAT )
  add_subdirectory( heat )
  list(APPEND TOOLBOXES "heat")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_HEAT AND FEELPP_TOOLBOXES_ENABLE_FLUIDMECHANICS )
  add_subdirectory( heatfluid )
  list(APPEND TOOLBOXES "heatfluid")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_ELECTRIC )
  add_subdirectory( electric )
  list(APPEND TOOLBOXES "electric")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_THERMOELECTRIC )
  add_subdirectory( thermoelectric )
  list(APPEND TOOLBOXES "thermoelectric")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MESHALE )
  add_subdirectory( alemesh )
  list(APPEND TOOLBOXES "alemesh")
endif()
if ( FEELPP_TOOLBOXES_ENABLE_MAXWELL )
  add_subdirectory( maxwell )
  list(APPEND TOOLBOXES "maxwell")
endif()

message(STATUS "[toolboxes] Toolboxes enabled : ${TOOLBOXES}")

# TOOLBOXES needs to be defined for pyfeelpp-toolboxes 
add_subdirectory( pyfeelpp-toolboxes )

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp-toolboxes.package.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/feelpp.version.cmake)
endif()


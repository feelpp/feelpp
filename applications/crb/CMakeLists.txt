###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-12
#
#  Copyright (C) 2009-2017 Université Joseph Fourier (Grenoble I)
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
project(crb)


# with this it is possible to compile applications/crb from an installed
# version of feel++. This is useful also for staged compilation for docker
if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  cmake_minimum_required(VERSION VERSION 3.1.3)
  find_package(Feel++ PATHS $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  message(STATUS "[feelpp] module path: ${CMAKE_MODULE_PATH}")
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
endif()

add_custom_target(crb)
set(FEELPP_CRB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(FEELPP_CRB_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})


include_directories(${FEELPP_CRB_SOURCE_DIR})
include_directories(${FEELPP_CRB_BINARY_DIR})

OPTION(FEELPP_ENABLE_CRB_ALL "enable all crb testcase" ON)

OPTION(FEELPP_ENABLE_CRB_EADS "enable EADS testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HEAT1D "enable HEAT1D testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HEAT2D "enable HEAT2D testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HEAT3D "enable HEAT3D testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_LINEARELASTICITY "enable LINEARELASTICITY testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_RBHEAT "enable RBHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_CABINEHEAT "enable CABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_FINALCABINEHEAT "enable FINALCABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_AD1 "enable AD1(Advection-Diffusion) testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_THERMALBLOCK "enable THERMALBLOCK(many parameters) testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HELMHOLTZ "enable Helmholtz testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HEATSINK "enable Heatsink testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_HEATSHIELD "enable Heatshield testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_BENCHMARKGREPL "enable BenchmarkGrepl testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_GEOHEATDIFFUSION "enable GeoHeatDiffusion testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_THERMALFIN "enable ThermalFin testcases" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_EEG "enable EEG testcase" OFF)
OPTION(FEELPP_ENABLE_CRB_CONVECTION "enable Convection testcase" ${FEELPP_ENABLE_CRB_ALL})
OPTION(FEELPP_ENABLE_CRB_STOKESDEIM "enable stokes with deim testcase" ${FEELPP_ENABLE_CRB_ALL})

set(FEELPP_CRB_LIBRARIES)
IF ( FEELPP_ENABLE_CRB_EADS )
  #add_subdirectory(eads)
ENDIF()

IF ( FEELPP_ENABLE_CRB_HEAT2D )
  add_subdirectory(heat2d)
ENDIF()

IF ( FEELPP_ENABLE_CRB_HEAT3D )
  add_subdirectory(heat3d)
ENDIF()

IF ( FEELPP_ENABLE_CRB_HEAT1D )
  add_subdirectory(heat1d)
  add_subdirectory(unsteady_heat1d)
ENDIF()

if ( FEELPP_ENABLE_CRB_LINEARELASTICITY )
    add_subdirectory(linearelasticity3d)
endif()
add_subdirectory(thermalbuilding)
add_subdirectory(opusheat)

IF ( FEELPP_ENABLE_CRB_RBHEAT )
  #add_subdirectory(Rbheat)
ENDIF()

IF ( FEELPP_ENABLE_CRB_CABINEHEAT )
  #add_subdirectory(Cabineheat)
ENDIF()

IF ( FEELPP_ENABLE_CRB_FINALCABINEHEAT )
  #add_subdirectory(Finalcabineheat)
ENDIF()

IF ( FEELPP_ENABLE_CRB_AD1 )
  add_subdirectory(advection-diffusion-2d-1)
ENDIF()


IF ( FEELPP_ENABLE_CRB_THERMALBLOCK )
  add_subdirectory(thermalblock)
ENDIF()


IF ( FEELPP_ENABLE_CRB_HELMHOLTZ )
add_subdirectory(helmholtz)
ENDIF()


if ( FEELPP_ENABLE_CRB_HEATSINK )
  add_subdirectory(heatsink-2d)
ENDIF()

if ( FEELPP_ENABLE_CRB_HEATSHIELD )
  add_subdirectory(heatshield)
ENDIF()

if ( FEELPP_ENABLE_CRB_BENCHMARKGREPL )
  add_subdirectory(BenchmarkGrepl)
ENDIF()

if ( FEELPP_ENABLE_CRB_GEOHEATDIFFUSION )
  add_subdirectory(GeoHeatDiffusion)
ENDIF()

if ( FEELPP_ENABLE_CRB_THERMALFIN )
  add_subdirectory(thermalfin)
ENDIF()

IF ( FEELPP_ENABLE_CRB_EEG )
  add_subdirectory(eeg)
ENDIF()

IF ( FEELPP_ENABLE_CRB_CONVECTION )
  add_subdirectory(convection)
ENDIF()

IF ( FEELPP_ENABLE_CRB_STOKESDEIM )
  add_subdirectory(stokesdeim)
ENDIF()



add_subdirectory(onlinerun)
add_subdirectory(minimal)
add_subdirectory(thermoelectric)
add_subdirectory(biotsavart)
add_subdirectory(grepldeim)


include(${FEELPP_CRB_SOURCE_DIR}/cmake/feelpp.crb.cmake)

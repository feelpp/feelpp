###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-12
#
#  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)
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
add_custom_target(crb)
set(FEELPP_CRB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(FEELPP_CRB_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
include(UseCRB)

OPTION(FEELPP_ENABLE_CRB_EADS "enable EADS testcase" ON)
OPTION(FEELPP_ENABLE_CRB_HEAT1D "enable HEAT1D testcase" ON)
OPTION(FEELPP_ENABLE_CRB_HEAT2D "enable HEAT2D testcase" ON)
OPTION(FEELPP_ENABLE_CRB_RBHEAT "enable RBHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_CABINEHEAT "enable CABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_FINALCABINEHEAT "enable FINALCABINEHEAT testcases" OFF)
OPTION(FEELPP_ENABLE_CRB_AD1 "enable AD1(Advection-Diffusion) testcase" ON)
OPTION(FEELPP_ENABLE_CRB_THERMALBLOCK "enable THERMALBLOCK(many parameters) testcase" ON)
OPTION(FEELPP_ENABLE_CRB_HELMHOLTZ "enable Helmholtz testcases" ON)
OPTION(FEELPP_ENABLE_CRB_HEATSINK "enable Heatsink testcases" ON)
OPTION(FEELPP_ENABLE_CRB_HEATSHIELD "enable Heatshield testcases" ON)
OPTION(FEELPP_ENABLE_CRB_BENCHMARKGREPL "enable BenchmarkGrepl testcases" ON)
OPTION(FEELPP_ENABLE_CRB_GEOHEATDIFFUSION "enable GeoHeatDiffusion testcases" ON)
OPTION(FEELPP_ENABLE_CRB_EEG "enable EEG testcase" OFF)

OPTION(FEELPP_ENABLE_CRB_ALL "enable all crb testcase" ON)
if (FEELPP_ENABLE_CRB_ALL)
  set(FEELPP_ENABLE_CRB_EADS ON)
  set(FEELPP_ENABLE_CRB_HEAT1D ON)
  set(FEELPP_ENABLE_CRB_HEAT2D ON)
  set(FEELPP_ENABLE_CRB_RBHEAT OFF)
  set(FEELPP_ENABLE_CRB_CABINEHEAT OFF)
  set(FEELPP_ENABLE_CRB_FINALCABINEHEAT OFF)
  set(FEELPP_ENABLE_CRB_HEATSINK ON)
  set(FEELPP_ENABLE_CRB_AD1 ON)
  set(FEELPP_ENABLE_CRB_THERMALBLOCK ON)
  set(FEELPP_ENABLE_CRB_HELMHOLTZ ON)
  set(FEELPP_ENABLE_CRB_HEATSHIELD ON)
  set(FEELPP_ENABLE_CRB_BENCHMARKGREPL ON)
  set(FEELPP_ENABLE_CRB_GEOHEATDIFFUSION ON)
  set(FEELPP_ENABLE_CRB_EEG OFF)
endif()


IF ( FEELPP_ENABLE_CRB_EADS )
  add_subdirectory(eads)
ENDIF()

IF ( FEELPP_ENABLE_CRB_HEAT2D )
  add_subdirectory(heat2d)
ENDIF()

IF ( FEELPP_ENABLE_CRB_HEAT1D )
  add_subdirectory(heat1d)
  add_subdirectory(unsteady_heat1d)
ENDIF()

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

IF ( FEELPP_ENABLE_CRB_EEG )
  add_subdirectory(eeg)
ENDIF()




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
add_custom_target(opus)
set(FEELPP_OPUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(FEELPP_OPUS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
include(UseCRB)


OPTION(FEELPP_OPUS_ENABLE_EADS "enable EADS testcase" OFF)
IF ( FEELPP_OPUS_ENABLE_EADS )
  add_subdirectory(eads)
ENDIF()

OPTION(FEELPP_OPUS_ENABLE_HEAT1D "enable HEAT1D testcase" OFF)
IF ( FEELPP_OPUS_ENABLE_HEAT1D )
  add_subdirectory(heat1d)
  add_subdirectory(unsteady_heat1d)
ENDIF()

OPTION(FEELPP_OPUS_ENABLE_AD1 "enable AD1(Advection-Diffusion) testcase" OFF)
IF ( FEELPP_OPUS_ENABLE_AD1 )
  add_subdirectory(advection-diffusion-2d-1)
ENDIF()

OPTION(FEELPP_OPUS_ENABLE_THERMALBLOCK "enable THERMALBLOCK(many parameters) testcase" OFF)
IF ( FEELPP_OPUS_ENABLE_THERMALBLOCK )
  add_subdirectory(thermalblock)
ENDIF()

OPTION(FEELPP_OPUS_ENABLE_HELMHOLTZ "enable Helmholtz testcases" OFF)
IF ( FEELPP_OPUS_ENABLE_HELMHOLTZ )
add_subdirectory(helmholtz)
ENDIF()

OPTION(FEELPP_OPUS_ENABLE_HEATSINK "enable Heatsink testcases" OFF)
if ( FEELPP_OPUS_ENABLE_HEATSINK )
  add_subdirectory(heatsink-2d)
endif()


###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-06-04
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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

#crb_add_executable(thermalblockpfem
#  thermalblockmain.cpp thermalblock.cpp  thermalblockgeo.cpp
#  LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG thermalblock.cfg )
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#crb_add_executable(thermalblockpfem thermalblock.cpp thermalblockpfem.cpp thermalblockgeo.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG thermalblockpfem.cfg )
#crb_add_executable(thermalblockcrb  thermalblock.cpp thermalblockcrb.cpp  thermalblockgeo.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG thermalblockcrb.cfg )

crb_add_model(thermalblock ThermalBlock HDRS thermalblock.hpp  LINK_LIBRARIES ${FEELPP_LIBRARIES}
  CFG thermalblock.cfg
  XML thermalblock.xml
  SCRIPTS thermalblock.m
  )
crb_add_model(thermalblockfree ThermalBlockFree HDRS thermalblock-operatorsfree.hpp  LINK_LIBRARIES ${FEELPP_LIBRARIES}
  CFG thermalblockfree.cfg
  XML thermalblockfree.xml
  SCRIPTS thermalblockfree.m
  )
crb_add_model(
  thermalblockminimalversion ThermalBlockMinimalVersion HDRS thermalblock-minimalversion.hpp  LINK_LIBRARIES ${FEELPP_LIBRARIES} 
  CFG thermalblockminimalversion.cfg 
  GEO thermalblock33.geo
  XML thermalblockminimalversion.xml
  SCRIPTS thermalblockminimalversion.m
  )

###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-11
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
feel_add_example( thermalfin thermalfin.cpp fin.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES})
feel_add_example( resistance resistance.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES})
feel_add_example( heatsink heatsink.cpp fin_sink.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heatsink.cfg)
feel_add_example( MCS_heat MCS_heat.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES})
feel_add_example( MCS_heat_time MCS_heat_time.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES})
#feel_add_example( Cabine_heat Cabine_heat.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES})
#feel_add_example( ground ground.cpp  LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG ground.cfg ground-realistic.cfg ) #ground1.cfg ground2.cfg ground3.cfg )
#feel_geo(${CMAKE_CURRENT_SOURCE_DIR}/ground.geo)
feel_add_example( Cabine Cabine_heat.cpp geometrie_cabine.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )

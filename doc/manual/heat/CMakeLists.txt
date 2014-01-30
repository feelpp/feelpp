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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application( thermalfin SRCS thermalfin.cpp LABEL doc GEO thermalfin.geo CFG thermalfin.cfg )
feelpp_add_application( resistance SRCS resistance.cpp LABEL doc)
feelpp_add_application( heatsink SRCS heatsink.cpp LABEL doc CFG heatsink.cfg GEO heatsink2d.geo)
feelpp_add_application( ground SRCS ground.cpp  LABEL doc CFG ground.cfg ground-realistic.cfg GEO ground.geo  ) #ground1.cfg ground2.cfg ground3.cfg )

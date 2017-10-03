###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-06-25
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
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


include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/feel/feelmodels/hifimagnet )

feelpp_add_application( thermoelectric_model_3D_V1T1_N1
                        SRCS thermoelectric_main.cpp ../feel/feelmodels/hifimagnet/options.cpp
                        LINK_LIBRARIES ${FEELPP_LIBRARIES}
                        DEFS FEELPP_DIM=3 FEELPP_ORDERV=1 FEELPP_ORDERT=1 FEELPP_GEO_ORDER=1
                        LABELS hifimagnet3D_P1_N1
                        INSTALL
                        )

feelpp_add_application( magnetostatic_mod_3D_P1_N1
                        SRCS magnetostatic_main2.cpp ../feel/feelmodels/hifimagnet/options.cpp
                        LINK_LIBRARIES ${FEELPP_LIBRARIES}
                        DEFS FEELPP_DIM=3 FEELPP_ORDER=1 FEELPP_G_ORDER=1
                        LABELS hifimagnet
                        NO_TEST
                        )

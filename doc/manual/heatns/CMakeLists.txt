###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-12-01
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
# Natural convection problem
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_library(
    lib_convection
    SHARED
    convection_other.cpp
    convection_run.cpp
    convection_residual.cpp
    convection_jacobian.cpp
    convection_jacobian1.cpp
    convection_jacobian2.cpp
    convection_lo.cpp
    convection_lo2.cpp )
target_link_libraries( lib_convection ${FEELPP_LIBRARIES}  )
add_executable( heatns_convection
                convection_main.cpp )
target_link_libraries( heatns_convection ${FEELPP_LIBRARIES} lib_convection )

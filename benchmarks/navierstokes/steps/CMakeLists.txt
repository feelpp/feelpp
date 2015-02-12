###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 08 Jan 2015
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


feelpp_add_application( steady_ns_2d SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_2d.cfg steady_ns_2d.bc DEFS FEELPP_DIM=2 FEELPP_ORDER_P=1)
feelpp_add_application( steady_ns_3d SRCS steady_ns.cpp INCLUDE_IN_ALL CFG steady_ns_3d.cfg steady_ns_3d.bc DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1)
install(TARGETS feelpp_steady_ns_2d DESTINATION bin/ COMPONENT Bin )
install(TARGETS feelpp_steady_ns_3d DESTINATION bin/ COMPONENT Bin )


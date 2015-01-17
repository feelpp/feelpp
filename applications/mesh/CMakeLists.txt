###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-06-11
#
#  Copyright (C) 2014 Feel++ Consortium
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

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR} )
feelpp_add_application( mesh3d SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 CONVEX=Simplex)
feelpp_add_application( mesh2d SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 CONVEX=Simplex)
feelpp_add_application( mesh1d SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=1 CONVEX=Simplex)
feelpp_add_application( mesh3d_hypercube SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 CONVEX=Hypercube)
feelpp_add_application( mesh2d_hypercube SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2 CONVEX=Hypercube)
feelpp_add_application( mesh1d_hypercube SRCS mesh.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=1 CONVEX=Hypercube)

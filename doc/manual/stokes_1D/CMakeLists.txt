###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-02-03
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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

feelpp_add_application(doc_stokes1D_P2P1 SRCS stokes1D.cpp LABELS doc GEO tige1D.geo cylinder3D.geo cylinder2D.geo INCLUDE_IN_ALL CFG stokes1D.cfg DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1)

feelpp_add_application(doc_stokes1D_P2P1_2D SRCS stokes1D.cpp LABELS doc GEO tige1D.geo cylinder3D.geo cylinder2D.geo INCLUDE_IN_ALL CFG stokes1D.cfg DEFS FEELPP_DIM=2 FEELPP_ORDER_P=1)

#feelpp_add_application(doc_mesh1D_P2P1 SRCS mesh1D.cpp LABELS doc GEO tige1D.geo INCLUDE_IN_ALL DEFS FEELPP_DIM=3 FEELPP_ORDER_P=1)

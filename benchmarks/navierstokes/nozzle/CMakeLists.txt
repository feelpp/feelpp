###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2012-12-16
#
#  Copyright (C) 2012 Feel++ Consortium
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

feelpp_add_application( ns_oseen_3d SRCS ns_oseen.cpp INCLUDE_IN_ALL GEO FDA2.geo CFG ns_oseen.cfg DEFS FEELPP_DIM=3 )

feelpp_add_application( ns_oseen_2d SRCS ns_oseen.cpp INCLUDE_IN_ALL
  GEO 2D_FDA_L.geo marche.geo cyl.geo step-erturk.geo
  CFG fda-2d.cfg marche.cfg step-erturk.cfg
  DEFS FEELPP_DIM=2 )

###  TEMPLATE.txt.tpl; coding: utf-8 ---
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2014-01-19
#
#  Copyright (C) 2011-2014 Feel++ Consortium
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
crb_add_model(geoheatdiffusion GeoHeatDiffusion HDRS geoheatdiffusion.hpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG geoheatdiffusion.cfg XML geoheatdiffusion.xml )

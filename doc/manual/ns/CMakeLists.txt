###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-09-30
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

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
#feelpp_add_application( backwardstep SRCS bs.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} DEFS FEELPP_NS_DIM=2)
feelpp_add_application( nsproj SRCS nsproj.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG nsproj.cfg GEO tubecylinder.geo)

add_subdirectory( drivencavity )

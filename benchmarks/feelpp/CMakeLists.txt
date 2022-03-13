###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 07 Apr 2019
#
#  Copyright (C) 2019 Feel++ Consortium
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
add_custom_target( benchs )

feelpp_add_application( env SRCS 01-env.cpp LINK_LIBRARIES  benchmark EXEC B_ENV )
feelpp_add_application( mesh SRCS 02-mesh.cpp LINK_LIBRARIES  benchmark EXEC B_MESH )
feelpp_add_application( space SRCS 03-space.cpp LINK_LIBRARIES  benchmark EXEC B_SPACE )
feelpp_add_application( integrate SRCS 04-integrate.cpp LINK_LIBRARIES  benchmark EXEC B_INTEGRATE )
feelpp_add_application( gmc SRCS gmc.cpp LINK_LIBRARIES  benchmark EXEC B_GMC )
feelpp_add_application( taskflow SRCS 99-taskflow.cpp LINK_LIBRARIES feelpp_taskflow )
add_dependencies( benchs ${B_ENV} ${B_MESH} ${B_SPACE} ${B_INTEGRATE} ${B_GMC})

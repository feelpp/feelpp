###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2011-04-23
#
#  Copyright (C) 2011 Universit� Joseph Fourier (Grenoble I)
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

mor_add_library( heat2d SRCS heat2d.cpp EXEC CRB_LIB_HEAT2D PLUGIN)
target_compile_definitions(${CRB_LIB_HEAT2D} PRIVATE FEELPP_MOR_PLUGIN_NAME=heat2d FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_HEAT2D})
install(FILES heat2d.hpp DESTINATION include/feelpp/mor/heat2d COMPONENT Devel)
set(FEELPP_CRB_HEAT2D_LIBRARIES ${CRB_LIB_HEAT2D} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_HEAT2D_LIBRARIES PARENT_SCOPE)

mor_add_model(heat2d Heat2D HDRS heat2d.hpp LINK_LIBRARIES ${CRB_LIB_HEAT2D} )

feelpp_add_testcase(heat2d PREFIX feelpp_crb CATEGORY CRB)
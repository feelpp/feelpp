###  TEMPLATE.txt.tpl; coding: utf-8 ---
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
mor_add_library( geoheatdiffusion SRCS geoheatdiffusion.cpp EXEC CRB_LIB_GEOHEAT PLUGIN)
target_compile_definitions(${CRB_LIB_GEOHEAT} PRIVATE FEELPP_MOR_PLUGIN_NAME=geoheatdiffusion FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_GEOHEAT})
install(FILES geoheatdiffusion.hpp DESTINATION include/feelpp/crb/geoheatdiffusion COMPONENT Devel)
set(FEELPP_CRB_GEOHEAT_LIBRARIES ${CRB_LIB_GEOHEAT} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_GEOHEAT_LIBRARIES PARENT_SCOPE)

mor_add_model(geoheatdiffusion GeoHeatDiffusion HDRS geoheatdiffusion.hpp LINK_LIBRARIES ${CRB_LIB_GEOHEAT})

feelpp_add_testcase(geoheatdiffusion PREFIX feelpp_crb CATEGORY CRB)
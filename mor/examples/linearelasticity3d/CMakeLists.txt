###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 12 Jun 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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


mor_add_library( linearelasticity3d SRCS linearelasticity3d.cpp EXEC CRB_LIB_LINEARELASTICITY3D)
install(FILES linearelasticity3d.hpp DESTINATION include/feelpp/crb/linearelasticity3d COMPONENT Devel)
set(FEELPP_CRB_LINEARELASTICITY3D_LIBRARIES ${CRB_LIB_LINEARELASTICITY3D} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_LINEARELASTICITY3D_LIBRARIES PARENT_SCOPE)

mor_add_model(linearelasticity3d LinearElasticity3d HDRS linearelasticity3d.hpp
  LINK_LIBRARIES ${CRB_LIB_LINEARELASTICITY3D} )

feelpp_add_testcase(linearelasticity3d PREFIX feelpp_crb CATEGORY CRB)

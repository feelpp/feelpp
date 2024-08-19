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

mor_add_library( heat1d SRCS heat1d.cpp EXEC CRB_LIB_HEAT1D PLUGIN EXPORT feelpp-mor-export-targets )
target_compile_definitions(${CRB_LIB_HEAT1D} PRIVATE FEELPP_MOR_PLUGIN_NAME=heat1d FEELPP_MOR_PLUGIN_LIBNAME=${CRB_LIB_HEAT1D})
install(FILES heat1d.hpp DESTINATION include/feelpp/mor/heat1d COMPONENT Devel)
set(FEELPP_CRB_HEAT1D_LIBRARIES ${CRB_LIB_HEAT1D} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_HEAT1D_LIBRARIES PARENT_SCOPE)


#mor_add_executable(heat1dpfem heat1dpfem.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heat1dpfem.cfg )
#mor_add_executable(heat1dpfem heat1dpfem.cpp LINK_LIBRARIES feel++_crb_crb CFG heat1dpfem.cfg )
#mor_add_executable(heat1dscm heat1dscm.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heat1dscm.cfg )
#mor_add_executable(heat1dcrb heat1dcrb.cpp LINK_LIBRARIES feel++_crb_crb ${FEELPP_LIBRARIES} CFG heat1dcrb.cfg )

#add_subdirectory(python)
#add_subdirectory(octave)

mor_add_model(heat1d Heat1D HDRS heat1d.hpp
  LINK_LIBRARIES ${CRB_LIB_HEAT1D}
  CFG heat1d.cfg XML heat1d.xml
  GEO heat1d.geo
  SCRIPTS heat1d.m)
#mor_add_python_module( heat1dcrb opusheat1dcrb_wrapper.cpp LINK_LIBRARIES feel++_crb_models SCRIPTS heat1d.py )

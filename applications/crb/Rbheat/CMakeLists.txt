###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-04-23
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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

#opus_add_executable(heat1dpfem heat1dpfem.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heat1dpfem.cfg )
#opus_add_executable(heat1dpfem heat1dpfem.cpp LINK_LIBRARIES feel++_opus_crb CFG heat1dpfem.cfg )
#opus_add_executable(heat1dscm heat1dscm.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG heat1dscm.cfg )
#opus_add_executable(heat1dcrb heat1dcrb.cpp LINK_LIBRARIES feel++_opus_crb ${FEELPP_LIBRARIES} CFG heat1dcrb.cfg )

#add_subdirectory(python)
#add_subdirectory(octave)


crb_add_model(Rbheat RbHeat HDRS Rbheat.hpp LINK_LIBRARIES ${FEELPP_LIBRARIES} CFG Rbheat.cfg XML Rbheat.xml SCRIPTS Rbheat.m)
#opus_add_python_module( heat1dcrb opusheat1dcrb_wrapper.cpp LINK_LIBRARIES feel++_opus_models SCRIPTS heat1d.py )

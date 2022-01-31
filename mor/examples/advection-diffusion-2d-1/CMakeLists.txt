###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2011-04-23
#
#  Copyright (C) 2011-2017 Universit� Joseph Fourier (Grenoble I)
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

mor_add_library( ad SRCS ad.cpp EXEC CRB_LIB_AD)
install(FILES ad.hpp DESTINATION include/feelpp/crb/ad COMPONENT Devel)
set(FEELPP_CRB_AD_LIBRARIES ${CRB_LIB_AD} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_AD_LIBRARIES PARENT_SCOPE)

#mor_add_executable(ad1pfem ad1pfem.cpp LINK_LIBRARIES feel++_crb_models ${FEELPP_LIBRARIES} CFG ad1pfem.cfg )
#mor_add_executable(ad1crb ad1crb.cpp LINK_LIBRARIES feel++_crb_crb ${FEELPP_LIBRARIES} CFG ad1crb.cfg)
mor_add_model(ad AdvectionDiffusion HDRS ad.hpp LINK_LIBRARIES ${CRB_LIB_AD} CFG ad.cfg XML ad.xml SCRIPTS ad.m )

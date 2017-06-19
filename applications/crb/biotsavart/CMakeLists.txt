###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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

include_directories(${CMAKE_SOURCE_DIR}/applications/crb/thermoelectric)

crb_add_library( biotsavart SRCS biotsavart.hpp biotsavart.cpp EXEC CRB_LIB_BIOTSAVART )
install( FILES biotsavart.hpp thermoelectric.hpp nl-thermoelectric.hpp
         DESTINATION include/feelpp/crb/biotsavart COMPONENT Devel)
set(FEELPP_CRB_BIOTSAVART_LIBRARIES ${CRB_LIB_BIOTSAVART} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_BIOTSAVART_LIBRARIES PARENT_SCOPE)

feelpp_add_application( biotsavartCRB
                        SRCS biotsavart_main.cpp biotsavart.hpp
                        LINK_LIBRARIES ${CRB_LIB_BIOTSAVART} )
feelpp_add_application( biotsavartCRBOpt
                        SRCS biotsavart_opt.cpp biotsavart.hpp
                        LINK_LIBRARIES ${CRB_LIB_BIOTSAVART} )

feelpp_add_testcase( bs-crb-quarterturn )

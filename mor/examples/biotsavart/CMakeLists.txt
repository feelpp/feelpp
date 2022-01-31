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

# mor_add_library( biotsavart SRCS biotsavart.hpp biotsavart.cpp LINK_LIBRARIES ${FEELPP_CRB_THERMOELECTRIC_LINEAR_LIBRARIES} EXEC CRB_LIB_BIOTSAVART )
# target_include_directories(${CRB_LIB_BIOTSAVART} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)
# set(FEELPP_CRB_BIOTSAVART_LIBRARIES ${CRB_LIB_BIOTSAVART} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_BIOTSAVART_LIBRARIES PARENT_SCOPE)

# mor_add_library( alphabiotsavart SRCS alphabiotsavart.hpp alphabiotsavart.cpp LINK_LIBRARIES ${FEELPP_CRB_ALPHATHERMOELECTRIC_LIBRARIES} EXEC CRB_LIB_ALPHABIOTSAVART )
# target_include_directories(${CRB_LIB_ALPHABIOTSAVART} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)
# set(FEELPP_CRB_ALPHABIOTSAVART_LIBRARIES ${CRB_LIB_ALPHABIOTSAVART} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ALPHABIOTSAVART_LIBRARIES PARENT_SCOPE)

mor_add_library( biotsavart-alpha-electric SRCS biotsavart-alpha-electric.hpp biotsavart-alpha-electric.cpp LINK_LIBRARIES ${FEELPP_CRB_ELECTRIC_ALPHA_LIBRARIES} EXEC CRB_LIB_BIOTSAVART_ALPHA_ELECTRIC )
target_include_directories(${CRB_LIB_BIOTSAVART_ALPHA_ELECTRIC} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)
set(FEELPP_CRB_BIOTSAVART_ALPHA_ELECTRIC_LIBRARIES ${CRB_LIB_BIOTSAVART_ALPHA_ELECTRIC} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_BIOTSAVART_ALPHA_ELECTRIC_LIBRARIES PARENT_SCOPE)

install( FILES biotsavart.hpp biotsavart-alpha-electric.hpp
         DESTINATION include/feelpp/crb/biotsavart COMPONENT Devel)

# feelpp_add_application( biotsavartCRB
#                         SRCS biotsavart_main.cpp biotsavart.hpp
#                         LINK_LIBRARIES ${CRB_LIB_BIOTSAVART}
#                         EXEC BIOTSAVART_APP
#                         )
# target_include_directories(${BIOTSAVART_APP} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)

# feelpp_add_application( biotsavartCRBOpt
#                         SRCS biotsavart_opt.cpp biotsavart.hpp
#                         LINK_LIBRARIES ${CRB_LIB_BIOTSAVART}
#                         EXEC BIOTSAVARTOPT_APP )
# target_include_directories(${BIOTSAVARTOPT_APP} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)

# feelpp_add_application( biotsavartCRBAlpha
#                         SRCS alphabiotsavart_main.cpp alphabiotsavart.hpp
#                         LINK_LIBRARIES ${CRB_LIB_ALPHABIOTSAVART}
#                         EXEC ALPHABIOTSAVART_APP
#                         )
# target_include_directories(${ALPHABIOTSAVART_APP} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)

feelpp_add_application( biotsavart-alpha-electric-app
                        SRCS biotsavart-alpha-electric_main.cpp biotsavart-alpha-electric.hpp
                        LINK_LIBRARIES ${CRB_LIB_BIOTSAVART_ALPHA_ELECTRIC}
                        EXEC BIOTSAVART_ALPHA_ELECTRIC_APP
                        )
target_include_directories(${BIOTSAVART_ALPHA_ELECTRIC_APP} PRIVATE ${PROJECT_SOURCE_DIR}/thermoelectric)


feelpp_add_testcase( bs-crb-quarterturn )
feelpp_add_testcase( bs-bancmesure )
feelpp_add_testcase( bs-alpha-hel )
feelpp_add_testcase( crb-bs-alpha-electric )


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

crb_add_library( electric SRCS electric.cpp EXEC CRB_LIB_ELECTRIC)
set(FEELPP_CRB_ELECTRIC_LIBRARIES ${CRB_LIB_ELECTRIC} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ELECTRIC_LIBRARIES PARENT_SCOPE)

crb_add_library( thermic SRCS thermic.cpp EXEC CRB_LIB_THERMIC)
set(FEELPP_CRB_THERMIC_LIBRARIES ${CRB_LIB_THERMIC} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMIC_LIBRARIES PARENT_SCOPE)

crb_add_library( thermoelectric-simple SRCS thermoelectric-simple.cpp EXEC CRB_LIB_THERMOELECTRIC_SIMPLE)
set(FEELPP_CRB_THERMOELECTRIC_SIMPLE_LIBRARIES ${CRB_LIB_THERMOELECTRIC_SIMPLE} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMOELECTRIC_SIMPLE_LIBRARIES PARENT_SCOPE)

# crb_add_library( thermoelectric SRCS thermoelectric.cpp EXEC CRB_LIB_THERMOELECTRIC)
# set(FEELPP_CRB_THERMOELECTRIC_LIBRARIES ${CRB_LIB_THERMOELECTRIC} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMOELECTRIC_LIBRARIES PARENT_SCOPE)

# crb_add_library( nlthermoelectric SRCS nlthermoelectric.cpp EXEC CRB_LIB_NLTHERMOELECTRIC)
# set(FEELPP_CRB_NLTHERMOELECTRIC_LIBRARIES ${CRB_LIB_NLTHERMOELECTRIC} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_NLTHERMOELECTRIC_LIBRARIES PARENT_SCOPE)

# need DEIM for composite space
# crb_add_library( alphathermoelectric SRCS alphathermoelectric.cpp EXEC CRB_LIB_ALPHATHERMOELECTRIC)
# set(FEELPP_CRB_ALPHATHERMOELECTRIC_LIBRARIES ${CRB_LIB_ALPHATHERMOELECTRIC} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ALPHATHERMOELECTRIC_LIBRARIES PARENT_SCOPE)

crb_add_library( electric-alpha SRCS electric-alpha.cpp EXEC CRB_LIB_ELECTRIC_ALPHA)
set(FEELPP_CRB_ELECTRIC_ALPHA_LIBRARIES ${CRB_LIB_ELECTRIC_ALPHA} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ELECTRIC_ALPHA_LIBRARIES PARENT_SCOPE)


# feelpp_add_application( alpha-thermoelectricCRB
#                         SRCS alphathermoelectric.hpp alphathermoelectric_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_ALPHATHERMOELECTRIC} INSTALL
#                         )

install( FILES electric.hpp thermic.hpp thermoelectric-simple.hpp electric-alpha.hpp # thermoelectric.hpp nlthermoelectric.hpp
         DESTINATION include/feelpp/crb/thermoelectric COMPONENT Devel)


feelpp_add_application( electric-app
                        SRCS electric.hpp electric_main.cpp
                        LINK_LIBRARIES ${CRB_LIB_ELECTRIC} INSTALL )
feelpp_add_application( thermic-app
                        SRCS thermic.hpp thermic_main.cpp
                        LINK_LIBRARIES ${CRB_LIB_THERMIC} INSTALL )
feelpp_add_application( thermoelectric-simple-app
                        SRCS thermoelectric-simple.hpp thermoelectric-simple_main.cpp
                        LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_SIMPLE} INSTALL )
# feelpp_add_application( thermoelectricCRB
#                         SRCS thermoelectric.hpp thermoelectric_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC} INSTALL )
# feelpp_add_application( nl-thermoelectricCRB
#                         SRCS nlthermoelectric.hpp nlthermoelectric_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_NLTHERMOELECTRIC} INSTALL )
feelpp_add_application( electric-alpha-app
                        SRCS electric-alpha_main.cpp electric-alpha.hpp
                        LINK_LIBRARIES ${CRB_LIB_ELECTRIC_ALPHA}
                        EXEC ELECTRIC_ALPHA_APP
                        )

feelpp_add_application( alphaElectricFE SRCS alphaelectricFE_main.cpp )
# feelpp_add_application( testthCRB SRCS test.cpp )
# feelpp_add_application( testthCRB SRCS test.cpp )
# feelpp_add_application( thermoelectric_geo_FE_2D
#                         SRCS thermoelectric_geo_FE.cpp DEFS FEELPP_DIM=2 )
# feelpp_add_application( thermoelectric_geo_FE_3D
#                         SRCS thermoelectric_geo_FE.cpp DEFS FEELPP_DIM=3 )

feelpp_add_testcase( th-crb-quarterturn CATEGORY CRB)
feelpp_add_testcase( te-bancmesure CATEGORY CRB)
feelpp_add_testcase( crb-electric-alpha CATEGORY CRB)
feelpp_add_testcase( crb-electric CATEGORY CRB)
feelpp_add_testcase( crb-thermic CATEGORY CRB)
feelpp_add_testcase( crb-thermoelectric-simple CATEGORY CRB)

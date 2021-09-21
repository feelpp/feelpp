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

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  find_package(Feel++ COMPONENTS Toolboxes PATHS ${FEELPP_DIR}/share/feelpp/feel/cmake/modules $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules )
  if(NOT FEELPP_FOUND)
    message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.")
  endif()
  if (NOT FEELPP_TOOLBOXES_FOUND)
    message(FATAL_ERROR "Feel++ was not installed with toolboxes")
  endif()
endif()

# how to FEELPP_INCLUDE_DIR
if (EXISTS ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
  MESSAGE(STATUS "Loading toolboxes headers from ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes")
  include_directories(${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
else()
  MESSAGE(STATUS "Loading toolboxes headers from ${CMAKE_SOURCE_DIR}/toolboxes and ${CMAKE_BINARY_DIR}/toolboxes")
  include_directories(${CMAKE_SOURCE_DIR}/toolboxes)
  include_directories(${CMAKE_BINARY_DIR}/toolboxes)
endif()

# crb_add_library( electric SRCS electric.cpp EXEC CRB_LIB_ELECTRIC)
# set(FEELPP_CRB_ELECTRIC_LIBRARIES ${CRB_LIB_ELECTRIC} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ELECTRIC_LIBRARIES PARENT_SCOPE)

# crb_add_library( thermic SRCS thermic.cpp LINK_LIBRARIES ${CRB_LIB_ELECTRIC} EXEC CRB_LIB_THERMIC)
# set(FEELPP_CRB_THERMIC_LIBRARIES ${CRB_LIB_THERMIC} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMIC_LIBRARIES PARENT_SCOPE)

# crb_add_library( thermoelectric-simple SRCS thermoelectric-simple.cpp EXEC CRB_LIB_THERMOELECTRIC_SIMPLE)
# set(FEELPP_CRB_THERMOELECTRIC_SIMPLE_LIBRARIES ${CRB_LIB_THERMOELECTRIC_SIMPLE} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMOELECTRIC_SIMPLE_LIBRARIES PARENT_SCOPE)

# crb_add_library( thermoelectric-linear SRCS thermoelectric-linear.cpp EXEC CRB_LIB_THERMOELECTRIC_LINEAR)
# set(FEELPP_CRB_THERMOELECTRIC_LINEAR_LIBRARIES ${CRB_LIB_THERMOELECTRIC_LINEAR} PARENT_SCOPE)
# set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMOELECTRIC_LINEAR_LIBRARIES PARENT_SCOPE)

crb_add_library( thermoelectric-nonlinear SRCS thermoelectric-nl_create.cpp thermoelectric-nl_assemble.cpp thermoelectric-nl_others.cpp EXEC CRB_LIB_THERMOELECTRIC_NONLINEAR)
set(FEELPP_CRB_THERMOELECTRIC_NONLINEAR_LIBRARIES ${CRB_LIB_THERMOELECTRIC_NONLINEAR} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_THERMOELECTRIC_NONLINEAR_LIBRARIES PARENT_SCOPE)

crb_add_library( electric-alpha SRCS electric-alpha.cpp EXEC CRB_LIB_ELECTRIC_ALPHA)
set(FEELPP_CRB_ELECTRIC_ALPHA_LIBRARIES ${CRB_LIB_ELECTRIC_ALPHA} PARENT_SCOPE)
set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_ELECTRIC_ALPHA_LIBRARIES PARENT_SCOPE)

if (FEELPP_HAS_GLPK)
  crb_add_library( biotsavart-electric-alpha SRCS biotsavart-electric-alpha.hpp biotsavart-electric-alpha.cpp LINK_LIBRARIES ${CRB_LIB_ELECTRIC_ALPHA} ${GLPK_LIBRARIES} EXEC CRB_LIB_BIOTSAVART_ELECTRIC_ALPHA )
  target_include_directories(feelpp_crb_biotsavart-electric-alpha PUBLIC ${GLPK_INCLUDE_DIR} )
  set(FEELPP_CRB_BIOTSAVART_ELECTRIC_ALPHA_LIBRARIES ${CRB_LIB_BIOTSAVART_ELECTRIC_ALPHA} PARENT_SCOPE)
  set(FEELPP_CRB_LIBRARIES ${FEELPP_CRB_LIBRARIES} FEELPP_CRB_BIOTSAVART_ELECTRIC_ALPHA_LIBRARIES PARENT_SCOPE)
endif()

#install( FILES electric.hpp thermic.hpp thermoelectric-simple.hpp thermoelectric-linear.hpp thermoelectric-nonlinear.hpp electric-alpha.hpp
install( FILES thermoelectric-nl.hpp electric-alpha.hpp biotsavart.hpp biotsavart-electric-alpha.hpp
         DESTINATION include/feelpp/crb/thermoelectric COMPONENT Devel)


# feelpp_add_application( electric-app
#                         SRCS electric.hpp electric_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_ELECTRIC} INSTALL )
# feelpp_add_application( thermic-app
#                         SRCS thermic.hpp thermic_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_THERMIC} INSTALL )
# feelpp_add_application( thermoelectric-simple-app
#                         SRCS thermoelectric-simple.hpp thermoelectric-simple_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_SIMPLE} INSTALL )
# feelpp_add_application( thermoelectric-linear-app
#                         SRCS thermoelectric-linear.hpp thermoelectric-linear_main.cpp
#                         LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_LINEAR} INSTALL )
# feelpp_add_application( saddle-thermoelectric-linear-app
#                         SRCS thermoelectric-linear.hpp thermoelectric-linear_main.cpp
#                         DEFS THERMOELECTRIC_SADDLEPOINT
#                         LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_LINEAR} INSTALL )
feelpp_add_application( thermoelectric-nonlinear-app
                        SRCS thermoelectric-nl.hpp thermoelectric-nl_main.cpp biotsavart-maxwell.hpp
                        LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_NONLINEAR} INSTALL )
feelpp_add_application( saddle-thermoelectric-nonlinear-app
                        SRCS thermoelectric-nl.hpp thermoelectric-nl_main.cpp
                        DEFS THERMOELECTRIC_SADDLEPOINT
                        LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_NONLINEAR} INSTALL )
feelpp_add_application( electric-alpha-app
                        SRCS electric-alpha_main.cpp electric-alpha.hpp
                        LINK_LIBRARIES ${CRB_LIB_ELECTRIC_ALPHA}
                        EXEC ELECTRIC_ALPHA_APP
                        )
feelpp_add_application( biotsavart-app
                        SRCS thermoelectric-nl.hpp biotsavart.hpp biotsavart_main.cpp
                        LINK_LIBRARIES ${CRB_LIB_THERMOELECTRIC_NONLINEAR} INSTALL )
if (FEELPP_HAS_GLPK)
  feelpp_add_application( biotsavart-electric-alpha-app
    SRCS electric-alpha.hpp biotsavart-electric-alpha.hpp biotsavart-electric-alpha_main.cpp
    LINK_LIBRARIES ${CRB_LIB_BIOTSAVART_ELECTRIC_ALPHA} Feelpp::feelpp_toolbox_electric_lib Feelpp::feelpp_modelcore )
endif()
# feelpp_add_application( test SRCS test.cpp )

# crb_add_model(thermoelectric-nonlinear ThermoElectric
# 	HDRS thermoelectric-nonlinear.hpp
# 	LINK_LIBRARIES ${FEELPP_LIBRARIES} ${CRB_LIB_THERMOELECTRIC_NONLINEAR}
#     )

# feelpp_add_testcase( crb-electric CATEGORY CRB)
# feelpp_add_testcase( crb-thermic CATEGORY CRB)
# feelpp_add_testcase( crb-thermoelectric-simple CATEGORY CRB)
# feelpp_add_testcase( crb-thermoelectric-linear CATEGORY CRB)
feelpp_add_testcase( crb-thermoelectric-nonlinear CATEGORY CRB)
feelpp_add_testcase( crb-electric-alpha CATEGORY CRB)
feelpp_add_testcase( crb-biotsavart CATEGORY CRB)
feelpp_add_testcase( bs-alpha-hel CATEGORY CRB)
feelpp_add_testcase( eq-test CATEGORY CRB)


if (FEELPP_HAS_GLPK)
  feelpp_add_application( test_eq SRCS test_eq.cpp LINK_LIBRARIES ${GLPK_LIBRARIES} )
  target_include_directories(feelpp_mor_test_eq PUBLIC ${GLPK_INCLUDE_DIR} )
endif()

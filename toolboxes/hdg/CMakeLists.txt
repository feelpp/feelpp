###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 30 Mar 2016
#
#  Copyright (C) 2016 Feel++ Consortium
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

set(FEELPP_HDG_PATH ${CMAKE_SOURCE_DIR}/toolboxes/feel/feelmodels/hdg )
#feelpp_add_application( hdg_darcy main.cpp )
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

if ( FEELPP_HAS_OMC AND FEELPP_HAS_FMILIB )
  add_subdirectory( coupledpoisson )
endif()

feelpp_add_application( hdg_poisson
  SRCS mixedpoisson_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib
  INSTALL MAN hdg_poisson
  TESTS )

feelpp_add_application( hdg_elasticity
  SRCS mixedelasticity_main.cpp
  INCLUDE_IN_ALL
  LINK_LIBRARIES feelpp_toolbox_hdg_elasticity_lib
  INSTALL MAN hdg_elasticity )

# feelpp_add_application( hdg_thermoelectric
#   SRCS thermoelectric.hpp thermoelectric_main.cpp
#   LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib
#   INSTALL
#   )

# feelpp_add_application( hdg_cvg_thermoelectric
#   SRCS thermoelectric.hpp cvg_thermoelectric.cpp
#   LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib
#   INSTALL
#   )

# feelpp_add_application( hdg_poroelastic
#   SRCS poroelastic_main.cpp
#   LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib feelpp_toolbox_hdg_elasticity_lib
#   INSTALL
#   )

#feelpp_add_application(
#  hdg_thermoelectricelasticity
#  SRCS thermoelectric-elasticity.hpp thermoelectric-elasticity_main.cpp
#  LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib feelpp_toolbox_hdg_elasticity_lib
#  )

if ( FEELPP_HAS_FMILIB )
  feelpp_add_application( hdg_coupledpoisson
    SRCS coupledmixedpoisson_main.cpp
    LINK_LIBRARIES feelpp_toolbox_hdg_poisson_lib
    INSTALL MAN hdg_coupledpoisson  )
endif()


foreach( _testcase IN ITEMS "poisson" "elasticity" "thermoelectric" "coupledpoisson" )
  feelpp_add_testcase(${_testcase} PREFIX feelpp_toolbox_hdg CATEGORY toolboxes/hdg)
endforeach()

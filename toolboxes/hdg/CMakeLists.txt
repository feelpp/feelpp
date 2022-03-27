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

if ( ( FEELPP_HAS_OMC OR FEELPP_HAS_MO2FMU ) AND FEELPP_HAS_FMI4CPP )
  add_subdirectory( cases/coupledpoisson )
endif()

feelpp_add_application( hdg_poisson
  SRCS mixedpoisson_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_lib
  INSTALL MAN hdg_poisson
  TESTS )

feelpp_add_application( hdg_elasticity
  SRCS mixedelasticity_main.cpp
  INCLUDE_IN_ALL
  LINK_LIBRARIES feelpp_toolbox_hdg_lib
  INSTALL MAN hdg_elasticity
  TESTS )

feelpp_add_application( hdg_thermoelectric
  SRCS thermoelectric_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_lib
  INSTALL TESTS
  )

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
feelpp_add_testcase(cases PREFIX feelpp_toolbox_hdg CATEGORY toolboxes/hdg)

if ( FEELPP_HAS_FMI4CPP )
  feelpp_add_application( hdg_coupledpoisson
    SRCS coupledmixedpoisson_main.cpp
    LINK_LIBRARIES feelpp_toolbox_hdg_lib
    INSTALL MAN hdg_coupledpoisson
     )

  message(STATUS "[toolbox] add tests feelpp_toolbox_hdg_coupledpoisson-linear")
  add_test( NAME feelpp_toolbox_hdg_coupledpoisson-linear-om COMMAND feelpp_toolbox_hdg_coupledpoisson --config-file ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupledpoisson/1-linear/linear.cfg
        --ts.time-final 0.25
        --hdg.poisson.json.patch "{ \"op\": \"replace\", \"path\": \"/BoundaryConditions/poisson/CouplingODEs/buffer/circuit\", \"value\": \"${CMAKE_CURRENT_BINARY_DIR}/cases/coupledpoisson/1-linear/test3d0d_linear/test3d0d_linear.fmu\" }" )

  # add_test( NAME feelpp_toolbox_hdg_coupledpoisson-linear-dy COMMAND feelpp_toolbox_hdg_coupledpoisson --config-file ${CMAKE_CURRENT_SOURCE_DIR}/cases/coupledpoisson/1-linear/linear.cfg --ts.time-final 0.25 --hdg.poisson.json-editions BoundaryConditions.flux.Coupling.buffer.circuit:${CMAKE_CURRENT_SOURCE_DIR}/cases/coupledpoisson/1-linear/test3d0d.fmu)
endif()




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

#feelpp_add_application( hdg_darcy main.cpp )
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

feelpp_add_application( mixed-poisson-model
  SRCS mixedpoisson_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_mixedpoisson_lib
  INSTALL )

feelpp_add_application( mixed-elasticity-model
  SRCS mixedelasticity_main.cpp
  INCLUDE_IN_ALL
  LINK_LIBRARIES feelpp_toolbox_hdg_mixedelasticity_lib
  INSTALL )

feelpp_add_application( thermoelectricHDG
  SRCS thermoelectric.hpp thermoelectric_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_mixedpoisson_lib
  INSTALL
  )

feelpp_add_application( poroelasticHDG
  SRCS poroelastic_main.cpp
  LINK_LIBRARIES feelpp_toolbox_hdg_mixedpoisson_lib feelpp_toolbox_hdg_mixedelasticity_lib
  INSTALL
  )
  endforeach()
endforeach()

feelpp_add_testcase(bitter CATEGORY hdg)
#feelpp_add_testcase(eye2brain_level1 CATEGORY hdg)
#feelpp_add_testcase(eye2brain_level2 CATEGORY hdg)
feelpp_add_testcase(cylinder CATEGORY hdg )
feelpp_add_testcase(testDarcy CATEGORY hdg )
feelpp_add_testcase(mixed-poisson CATEGORY hdg )
feelpp_add_testcase(mixed-elasticity CATEGORY hdg )
feelpp_add_testcase(mp-red CATEGORY hdg )
feelpp_add_testcase(HL-31)

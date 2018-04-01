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

foreach(D 3)
    foreach(P RANGE 0 ${FEELPP_INSTANTIATION_ORDER_MAX} ) #3
        feelpp_add_application( mixed-poisson-model_${D}DP${P}
                                SRCS mixedpoisson_main.cpp ${FEELPP_HDG_PATH}/options.cpp
                                LINK_LIBRARIES feelpp_modelcore
                                DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} FEELPP_MODELS_HDG_NOEXTERN
                                INSTALL )
  endforeach()
endforeach()

foreach(D 3)
  foreach(P  RANGE 0 ${FEELPP_INSTANTIATION_ORDER_MAX} ) #3
	foreach (N 1)
	  feelpp_add_application( mixed-elasticity-model_${D}DP${P}_G${N} 
							SRCS mixedelasticity_main.cpp ${FEELPP_HDG_PATH}/options.cpp
							DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} FEELPP_GEO_ORDER=${N} FEELPP_MODELS_HDG_NOEXTERN
							INCLUDE_IN_ALL
							LINK_LIBRARIES feelpp_modelcore
                            INSTALL )
	endforeach()
  endforeach()
endforeach()

foreach(D 2 3)
  foreach(P RANGE 0 ${FEELPP_INSTANTIATION_ORDER_MAX} ) #3
    feelpp_add_application( thermoelectricHDG_${D}DP${P}
                            SRCS thermoelectric.hpp thermoelectric_main.cpp ${FEELPP_HDG_PATH}/options.cpp
                            DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} FEELPP_MODELS_HDG_NOEXTERN
                            LINK_LIBRARIES feelpp_modelcore
                            INSTALL
                            )
  endforeach()
endforeach()

foreach(D 3)
  foreach(P 1 ) #3
    feelpp_add_application( poroelastic-model_${D}DP${P}
      SRCS poroelastic_main.cpp ${FEELPP_HDG_PATH}/options.cpp
      LINK_LIBRARIES feelpp_modelcore
      DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} FEELPP_MODELS_HDG_NOEXTERN
      INSTALL )
  endforeach()
endforeach()


feelpp_add_testcase(poisson PREFIX hdg CATEGORY HDG )
feelpp_add_testcase(elasticity PREFIX hdg CATEGORY HDG )
feelpp_add_testcase(poroelasticity PREFIX hdg CATEGORY HDG )


feelpp_add_testcase(bitter CATEGORY hdg)
#feelpp_add_testcase(eye2brain_level1 CATEGORY hdg)
#feelpp_add_testcase(eye2brain_level2 CATEGORY hdg)
feelpp_add_testcase(cylinder CATEGORY hdg )
feelpp_add_testcase(testDarcy CATEGORY hdg )
feelpp_add_testcase(mixed-poisson CATEGORY hdg )
feelpp_add_testcase(mixed-elasticity CATEGORY hdg )
feelpp_add_testcase(mp-red CATEGORY hdg )


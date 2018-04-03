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

foreach(D 3)
    foreach(P RANGE 0 4 0) #3 
        feelpp_add_application( mixed-poisson-model_${D}DP${P}
                                SRCS mixedpoisson_main.cpp
                                LINK_LIBRARIES feelpp_modelcore
                                DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P}
                                INSTALL )
  endforeach()
endforeach()

foreach(D 2 3)
  foreach(P  RANGE 0 4 0) #3 
	foreach (N RANGE 1 2 1)
	  feelpp_add_application( mixed-elasticity-model_${D}DP${P}_G${N} 
							SRCS mixedelasticity_main.cpp 
							DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P} FEELPP_GEO_ORDER=${N}
							INCLUDE_IN_ALL
							LINK_LIBRARIES feelpp_modelcore
                            INSTALL )
	endforeach()
  endforeach()
endforeach()

foreach(D 2 3)
  foreach(P RANGE 0 4 0) #3 
    feelpp_add_application( thermoelectricHDG_${D}DP${P}
                            SRCS thermoelectric.hpp thermoelectric_main.cpp
                            DEFS FEELPP_DIM=${D} FEELPP_ORDER=${P}
                            LINK_LIBRARIES feelpp_modelcore
                            INSTALL
                            )
  endforeach()
endforeach()


feelpp_add_testcase(bitter CATEGORY hdg)
feelpp_add_testcase(cylinder CATEGORY hdg )
feelpp_add_testcase(testDarcy CATEGORY hdg )
feelpp_add_testcase(mixed-poisson CATEGORY hdg )
feelpp_add_testcase(mixed-elasticity CATEGORY hdg )
feelpp_add_testcase(mp-red CATEGORY hdg )

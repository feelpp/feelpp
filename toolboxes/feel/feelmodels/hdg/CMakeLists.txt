###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 05 Sep 2018
#
#  Copyright (C) 2018 Feel++ Consortium
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

feelpp_toolboxes_add_library( hdgbase SRCS enums.cpp mixedpoissonboundaryconditions.cpp LINK_LIBRARIES feelpp_modelcore )

foreach(_DIM 2 3)
  math(EXPR _P_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDPOISSON ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 ${MAX_ORDER} )
    foreach(_G_ORDER RANGE 1 2 ) #${MAX_ORDER})
      genLibHdg(
        DIM     ${_DIM}
        P_ORDER ${_P_ORDER}
        GEO_ORDER ${_G_ORDER}
        POLYSET Vectorial
        )
      set(FEELPP_MODELS_LIBRARIES_HDG ${FEELPP_MODELS_LIBRARIES_HDG} ${HDG_LIB_NAME} )
    endforeach()
  endforeach()
endforeach()


foreach(_DIM 2 3)
  math(EXPR _P_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDELASTICITY ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 ${MAX_ORDER} )
    foreach(_G_ORDER RANGE 1 2) #${MAX_ORDER})
      genLibHdg(
        DIM     ${_DIM}
        P_ORDER ${_P_ORDER}
        GEO_ORDER ${_G_ORDER}
        POLYSET Tensor2Symm
        )
      set(FEELPP_MODELS_LIBRARIES_HDG ${FEELPP_MODELS_LIBRARIES_HDG} ${HDG_LIB_NAME} )
    endforeach()
  endforeach()
endforeach()

# TODO create interface lib
feelpp_toolboxes_add_library(hdg_lib SRCS libmixedpoisson.cpp LINK_LIBRARIES ${FEELPP_MODELS_LIBRARIES_HDG})

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/hdg COMPONENT Devel)

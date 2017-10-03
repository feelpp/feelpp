###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-06-25
#
#  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
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


include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

  
foreach(D 2 3)
  # geo order 1
  foreach(T  RANGE 1 ${PMAX} 1)
    foreach(V  RANGE ${T} ${PMAX} 1)
      if( (D MATCHES 3) AND ((T MATCHES 1) AND (V MATCHES 1)) )
        feelpp_add_application( thermoelectric_model_${D}D_V${V}T${T}_N1
          SRCS thermoelectric_main.cpp thermoelectric.hpp ${OPTIONS_SRC} ${THERMOELECTRIC_SRCS}
          CFG thermoelectric.cfg
          LINK_LIBRARIES ${FEELPP_LIBRARIES}
          DEFS FEELPP_DIM=${D} FEELPP_ORDERV=${V} FEELPP_ORDERT=${T} FEELPP_GEO_ORDER=1
          # LABELS hifimagnet${D}D_P${T}_N1
          INSTALL
          )
      else ()
        feelpp_add_application( thermoelectric_model_${D}D_V${V}T${T}_N1
          SRCS thermoelectric_main.cpp thermoelectric.hpp ${OPTIONS_SRC} ${THERMOELECTRIC_SRCS}
          CFG thermoelectric.cfg
          LINK_LIBRARIES ${FEELPP_LIBRARIES}
          DEFS FEELPP_DIM=${D} FEELPP_ORDERV=${V} FEELPP_ORDERT=${T} FEELPP_GEO_ORDER=1
          LABELS ${D}D_P${T}_N1
          )
      endif()
    endforeach()
  endforeach()

  # iso-parametric
  math(EXPR _D_ORDER '${PMAX}-1')
  feelpp_min( MAX_ORDER ${PMAX} ${_D_ORDER})
  foreach(T  RANGE 2 ${MAX_ORDER} 1)
    foreach(V  RANGE ${T} ${MAX_ORDER} 1)
        feelpp_add_application( thermoelectric_model_${D}D_V${V}T${T}_N${T}
          SRCS thermoelectric_main.cpp thermoelectric.hpp ${OPTIONS_SRC} ${THERMOELECTRIC_SRCS}
          CFG thermoelectric.cfg
          LINK_LIBRARIES ${FEELPP_LIBRARIES}
          DEFS FEELPP_DIM=${D} FEELPP_ORDERV=${V} FEELPP_ORDERT=${T} FEELPP_GEO_ORDER=${T}
          LABELS  ${D}D_P${T}_N${T}
          )
    endforeach()
  endforeach()

endforeach()

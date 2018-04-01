###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 18 Jan 2017
#
#  Copyright (C) 2017 Feel++ Consortium
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


# geo order 1 (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_max( MAX_ORDER 1 ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 1 ${MAX_ORDER} )
    genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

# isoparametric (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 2 ${MAX_ORDER} )
    genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      GEO_ORDER ${_POLY_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME} )
  endforeach()
endforeach()

# vectorial geo order 1 (non-periodic)
foreach(_DIM 2 3)
    math(EXPR _POLY_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}-1')
    feelpp_max( MAX_ORDER 1 ${_POLY_ORDER})
  message(STATUS "[feelpp] ADVECTION ${_DIM}D PNVecG1 N=1..${MAX_ORDER}")
  foreach( _POLY_ORDER RANGE 1 ${MAX_ORDER} )
      genLibAdvection(
      DIM     ${_DIM}
      POLY_ORDER ${_POLY_ORDER}
      POLY_SET Vectorial
      GEO_ORDER 1
      )
  set(FEELPP_MODELS_LIBRARIES_ADVECTION  ${FEELPP_MODELS_LIBRARIES_ADVECTION} ${LIBBASE_NAME}  )
  endforeach()
endforeach()

add_library(feelpp_model_advection SHARED libadvection.cpp )
target_link_libraries(feelpp_model_advection ${FEELPP_MODELS_LIBRARIES_ADVECTION})

install(TARGETS feelpp_model_advection DESTINATION lib/ COMPONENT Libs)
FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/advection COMPONENT Devel)

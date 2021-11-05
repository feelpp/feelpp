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

foreach(_DIM 2 3)
  math(EXPR _P_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDPOISSON ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 ${MAX_ORDER} )
    foreach( _G_ORDER RANGE 1 2 )
      genLibHdg(
        DIM     ${_DIM}
        P_ORDER ${_P_ORDER}
        GEO_ORDER ${_G_ORDER}
        POLYSET Scalar
        )
      set(FEELPP_MODELS_LIBRARIES_HDG ${FEELPP_MODELS_LIBRARIES_HDG} ${HDG_LIB_NAME} )
    endforeach()
  endforeach()
endforeach()
foreach(_DIM 2 3)
  math(EXPR _P_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDELASTICITY ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 1)#${MAX_ORDER} )
    foreach( _G_ORDER RANGE 1 1)#2 )
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

add_library(feelpp_toolbox_hdg_lib SHARED libmixedpoisson.cpp enums.cpp )
add_library(Feelpp::feelpp_toolbox_hdg_lib ALIAS feelpp_toolbox_hdg_lib )
set_target_properties(feelpp_toolbox_hdg_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_hdg_lib" )
target_link_libraries(feelpp_toolbox_hdg_lib ${FEELPP_MODELS_LIBRARIES_HDG})
set_property(TARGET feelpp_toolbox_hdg_lib PROPERTY MACOSX_RPATH ON)
  set_target_properties(feelpp_toolbox_hdg_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_hdg_lib EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/hdg COMPONENT Devel)

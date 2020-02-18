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
  math(EXPR _P_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDPOISSON ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 ${MAX_ORDER} )
    foreach( _G_ORDER RANGE 1 2 )
      genLibMixedPoisson(
        DIM     ${_DIM}
        P_ORDER ${_P_ORDER}
        GEO_ORDER ${_G_ORDER}
        )
      set(FEELPP_MODELS_LIBRARIES_MIXEDPOISSON ${FEELPP_MODELS_LIBRARIES_MIXEDPOISSON} ${MIXEDPOISSON_LIB_NAME} )
    endforeach()
  endforeach()
endforeach()

add_library(feelpp_toolbox_hdg_poisson_lib SHARED libmixedpoisson.cpp enums.cpp )
add_library(Feelpp::feelpp_toolbox_hdg_poisson_lib ALIAS feelpp_toolbox_hdg_poisson_lib )
set_target_properties(feelpp_toolbox_hdg_poisson_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_hdg_poisson_lib" )
target_link_libraries(feelpp_toolbox_hdg_poisson_lib ${FEELPP_MODELS_LIBRARIES_MIXEDPOISSON})
set_property(TARGET feelpp_toolbox_hdg_poisson_lib PROPERTY MACOSX_RPATH ON)
  set_target_properties(feelpp_toolbox_hdg_poisson_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_hdg_poisson_lib
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Libs EXPORT feelpp-toolboxes-targets)

foreach(_DIM 2 3)
  math(EXPR _P_ORDER '${FEELPP_INSTANTIATION_ORDER_MAX}')
  feelpp_max( MAX_ORDER 1 ${_P_ORDER})
  message(STATUS "[feelpp] MIXEDELASTICITY ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _P_ORDER RANGE 1 ${MAX_ORDER} )
    genLibMixedElasticity(
      DIM     ${_DIM}
      P_ORDER ${_P_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_MIXEDELASTICITY ${FEELPP_MODELS_LIBRARIES_MIXEDELASTICITY} ${MIXEDELASTICITY_LIB_NAME} )
  endforeach()
endforeach()

add_library(feelpp_toolbox_hdg_elasticity_lib SHARED libmixedelasticity.cpp )
add_library(Feelpp::feelpp_toolbox_hdg_elasticity_lib ALIAS feelpp_toolbox_hdg_elasticity_lib )
set_target_properties(feelpp_toolbox_hdg_elasticity_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_hdg_elasticity_lib" )
target_link_libraries(feelpp_toolbox_hdg_elasticity_lib ${FEELPP_MODELS_LIBRARIES_MIXEDELASTICITY})
set_property(TARGET feelpp_toolbox_hdg_elasticity_lib PROPERTY MACOSX_RPATH ON)
  set_target_properties(feelpp_toolbox_hdg_elasticity_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_hdg_elasticity_lib
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/hdg COMPONENT Devel)

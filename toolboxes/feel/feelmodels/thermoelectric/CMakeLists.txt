###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Vincent Chabannes
#       Date: 25 Nov 2017
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


# geo order 1
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_max( MAX_ORDER 1 ${_T_ORDER})
  message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNG1 N=1..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 1 ${MAX_ORDER} )
    genLibThermoelectric(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      P_ORDER ${_T_ORDER}
      GEO_ORDER 1
      )
    set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME}  )
  endforeach()
endforeach()

# isoparametric
foreach(_DIM 2 3)
  math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
  feelpp_min( MAX_ORDER ${FEELPP_MESH_MAX_ORDER} ${_T_ORDER})
  message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNGN N=2..${MAX_ORDER}")
  foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
    genLibThermoelectric(
      DIM     ${_DIM}
      T_ORDER ${_T_ORDER}
      P_ORDER ${_T_ORDER}
      GEO_ORDER ${_T_ORDER}
      )
    set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME} )
  endforeach()
endforeach()

if( ${FEELPP_MESH_MAX_ORDER} GREATER 1 )
  foreach(_DIM 2 3)
    math(EXPR _T_ORDER "${FEELPP_INSTANTIATION_ORDER_MAX}")
    feelpp_max( MAX_ORDER 1 ${_T_ORDER})
    message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D PNG2 N=2..${MAX_ORDER}")
    foreach( _T_ORDER RANGE 2 ${MAX_ORDER} )
      genLibThermoelectric(
        DIM     ${_DIM}
        T_ORDER ${_T_ORDER}
        P_ORDER ${_T_ORDER}
        GEO_ORDER 2
        )
      set(FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC} ${THERMOELECTRIC_LIB_NAME} )
    endforeach()
  endforeach()
endif()


add_library(feelpp_toolbox_thermoelectric_lib SHARED libthermoelectric.cpp )
add_library(Feelpp::feelpp_toolbox_thermoelectric_lib ALIAS feelpp_toolbox_thermoelectric_lib )  # to match exported target
set_target_properties(feelpp_toolbox_thermoelectric_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_thermoelectric_lib")
target_link_libraries(feelpp_toolbox_thermoelectric_lib ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC})
set_property(TARGET feelpp_toolbox_thermoelectric_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_thermoelectric_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_thermoelectric_lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Libs EXPORT feelpp-toolboxes-targets)

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/thermoelectric COMPONENT Devel)

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


# Keep only P1/G1 at compile time; all other orders are runtime values.
foreach(_DIM 2 3)
  foreach(_ORDER 1 Dynamic)
    message(STATUS "[feelpp] THERMOELECTRIC ${_DIM}D P${_ORDER}G${_ORDER}")
    genLibThermoelectric(
      DIM ${_DIM}
      T_ORDER ${_ORDER}
      P_ORDER ${_ORDER}
      GEO_ORDER ${_ORDER}
      )
    list(APPEND FEELPP_MODELS_LIBRARIES_THERMOELECTRIC ${THERMOELECTRIC_LIB_NAME})
  endforeach()
endforeach()

add_library(feelpp_toolbox_thermoelectric_lib SHARED libthermoelectric.cpp )
add_library(Feelpp::feelpp_toolbox_thermoelectric_lib ALIAS feelpp_toolbox_thermoelectric_lib )  # to match exported target
set_target_properties(feelpp_toolbox_thermoelectric_lib PROPERTIES OUTPUT_NAME "feelpp_toolbox_thermoelectric_lib")
target_link_libraries(feelpp_toolbox_thermoelectric_lib ${FEELPP_MODELS_LIBRARIES_THERMOELECTRIC})
set_property(TARGET feelpp_toolbox_thermoelectric_lib PROPERTY MACOSX_RPATH ON)
set_target_properties(feelpp_toolbox_thermoelectric_lib PROPERTIES VERSION ${FEELPP_TOOLBOXES_SHARED_VERSION} SOVERSION ${FEELPP_TOOLBOXES_SHARED_SOVERSION})
install(TARGETS feelpp_toolbox_thermoelectric_lib EXPORT feelpp-toolboxes-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

FILE(GLOB files "*.hpp")
INSTALL(FILES ${files} DESTINATION include/feelpp/toolboxes/feel/feelmodels/thermoelectric COMPONENT Devel)

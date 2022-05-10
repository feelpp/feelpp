###  CMakeLists.txt ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2009-11-13
#
#  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)
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
# set(CMAKE_CXX_FLAGS "-g -O0" )

set(SRC toolboxmor.cpp crbdb.cpp crbmodeldb.cpp pfemapp.cpp opusapp.cpp crbplugin_interface.cpp crbmodelproperties.cpp crbmodelparameters.cpp crbmodeloutputs.cpp )

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
add_library( feelpp_mor SHARED ${SRC} )
add_library(Feelpp::feelpp_mor ALIAS feelpp_mor )
set_target_properties(feelpp_mor PROPERTIES VERSION 1 SOVERSION 1)
target_link_libraries( feelpp_mor PUBLIC Feelpp::feelpp   )
target_include_directories(feelpp_mor PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/mor>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/mor>
  $<INSTALL_INTERFACE:include/feelpp> )  
  
#
# Install
#
INSTALL(FILES ${HDR} DESTINATION include/feelpp/feel/feelmor COMPONENT Devel)
INSTALL(TARGETS feelpp_mor EXPORT feelpp-mor-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )


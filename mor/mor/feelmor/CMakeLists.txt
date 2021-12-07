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
set(SRC crbdb.cpp crbmodeldb.cpp pfemapp.cpp opusapp.cpp crbplugin_interface.cpp toolboxmor.cpp options.cpp )


add_library( feelpp_mor SHARED ${SRC} )
add_library(Feelpp::feelpp_mor ALIAS feelpp_mor )
set_target_properties(feelpp_mor PROPERTIES VERSION 1 SOVERSION 1)
target_link_libraries( feelpp_mor PUBLIC Feelpp::feelpp   )
INSTALL(TARGETS feelpp_mor EXPORT feelpp-mor-export-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
target_include_directories(feelpp_mor PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
  $<INSTALL_INTERFACE:include/feelpp> )  
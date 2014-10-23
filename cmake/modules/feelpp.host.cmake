###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-17
#
#  Copyright (C) 2014 Feel++ Consortium
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
# OSX
if(APPLE)
  if ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
    OPTION( FEELPP_ENABLE_MACOSX_RPATH "Enables MACOSX Rpath feature" OFF )
    if  ( FEELPP_ENABLE_MACOSX_RPATH )
      set(CMAKE_MACOSX_RPATH ON)
      set(CMAKE_SKIP_BUILD_RPATH FALSE)
      set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
      set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
      set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
      message(STATUS "MACOSX RPATH enabled (CMAKE Version 3)")
    else()
      message(STATUS "MACOSX RPATH disabled (CMAKE Version 3)")
    endif()
  endif()
endif()

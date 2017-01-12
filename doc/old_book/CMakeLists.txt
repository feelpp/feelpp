###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-18
#
#  Copyright (C) 2013-2015 Feel++ Consortium
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
cmake_minimum_required(VERSION 2.8)
project(tumor)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  cmake_minimum_required(VERSION 2.8)
  FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH  /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
  if ( FEELPP_CMAKE_MODULES )
    set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
  else()
    message(FATAL_ERROR "Feel++ does not seem to have been installed on this platform")
  endif()
  Find_Package(Feel++)
endif()
# add_custom_target(check-tumor COMMAND "ctest" ARGS "-R" "tumor")

add_subdirectory( Tutorial/code )


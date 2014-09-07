###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-05-25
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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


# add_executable(adgenerate adgenerate.cpp)
# target_link_libraries(adgenerate qt-mt)


# add_custom_target(hppfiles WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ )
# add_dependencies(hppfiles adgenerate)
# add_custom_command(
#   TARGET hppfiles
#   POST_BUILD
#   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/adgenerate .)


add_executable(test_ad EXCLUDE_FROM_ALL test_ad.cpp )
target_link_libraries(test_ad feelpp ${FEELPP_LIBRARIES})


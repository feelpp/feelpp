# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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

add_subdirectory (feelcore)
add_subdirectory (feelalg)
add_subdirectory (feelmesh)
#add_subdirectory (feelpoly)
add_subdirectory (feeldiscr)
add_subdirectory (feelfilters)
add_subdirectory (feelmaterial)
add_subdirectory (feelvf )
add_subdirectory (feelsystem)

add_library(  feel SHARED options.cpp )
set_target_properties( feel PROPERTIES VERSION ${FEEL_SHARED_VERSION} SOVERSION ${FEEL_SHARED_SOVERSION})

target_link_libraries (feel feelmaterial feelvf feelfilters feeldiscr feelalg feelmesh feelcore )
set_target_properties(feel PROPERTIES VERSION ${FEEL_SHARED_VERSION} SOVERSION ${FEEL_SHARED_SOVERSION})

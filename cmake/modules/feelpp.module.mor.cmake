###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 23 Aug 2018
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

#
# mor
#

if ( FEELPP_ENABLE_MOR )
  if ( EXISTS ${CMAKE_SOURCE_DIR}/mor/CMakeLists.txt )
    SET(FEELPP_HAS_MOR 1)
    SET(FEELPP_ENABLED_MODULES "${FEELPP_ENABLED_MODULES} MOR" )
  else()
    MESSAGE(WARNING "[feelpp] Mor was not found on your system. Either install it or set FEELPP_ENABLE_MOR to OFF.")
  endif() 
endif()
if ( NOT FEELPP_HAS_MOR )
  SET(FEELPP_DISABLED_MODULES "${FEELPP_DISABLED_MODULES} MOR" )
endif()

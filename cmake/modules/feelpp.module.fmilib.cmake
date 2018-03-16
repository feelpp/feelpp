###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 02 Feb 2018
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

find_package( FMILib )

if( FMILIB_FOUND )
  message( STATUS "[fmilib] include dir : ${FMILIB_INCLUDE_DIR}")
  message( STATUS "[fmilib] library :  ${FMILIB_LIBRARIES}" )
  include_directories( ${FMILIB_INCLUDE_DIR} )
  set( FEELPP_HAS_FMILIB 1 )
  set(FEELPP_LIBRARIES ${FMILIB_LIBRARIES} ${FEELPP_LIBRARIES})
endif()



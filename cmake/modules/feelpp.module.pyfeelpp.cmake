###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 19 Jun 2017
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

#
# pyfeelpp
#

if ( FEELPP_ENABLE_PYFEELPP )
  if ( EXISTS ${CMAKE_SOURCE_DIR}/pyfeelpp/CMakeLists.txt )
    SET(FEELPP_HAS_PYFEELPP 1)
    SET(FEELPP_ENABLED_MODULES "${FEELPP_ENABLED_MODULES} PyFeelpp" )
    ADD_DEFINITIONS( -DFEELPP_HAS_PYFEELPP )
  else()
    MESSAGE(WARNING "[feelpp] PyFeelpp was not found on your system. Either install it or set FEELPP_ENABLE_PYFEELPP to OFF.")
  endif() 
endif()

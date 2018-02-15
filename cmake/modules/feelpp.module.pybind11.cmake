###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 10 Jul 2017
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
# PyBind11
#
option( FEELPP_ENABLE_PYBIND11 "Enable PYBIND11" ON )

if ( FEELPP_ENABLE_PYBIND11 )
  feelppContribPrepare( pybind11 )
  if( FEELPP_CONTRIB_PREPARE_SUCCEED )
    #add_subdirectory(pybind11)
    include_directories(${FEELPP_SOURCE_DIR}/contrib/pybind11/include)
    set(FEELPP_HAS_PYBIND11 1)
    set(FEELPP_ENABLED_OPTIONS "${FEELPP_ENABLED_OPTIONS} PyBind11" )
    #add_definitions( -DFEELPP_HAS_PYBIND11 )
  endif()
endif()

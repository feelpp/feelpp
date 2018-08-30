###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 30 Aug 2018
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

if(NOT FEELPP_ENABLE_BUILD_STATIC)
	set(CMAKE_C++_CREATE_SHARED_LIBRARY 1)
endif()

# feel++ lib
ADD_SUBDIRECTORY ( feel )

OPTION(FEELPP_ENABLE_QUICKSTART "enable Feel++ quickstart" ON)
OPTION(FEELPP_ENABLE_PYFEELPP_LIBFEELPP "enable PyFeel++" ON)

add_subdirectory(tools)
include(feelpp.module.quickstart)
if ( FEELPP_HAS_QUICKSTART )
  add_custom_target(quickstart)
  add_subdirectory(quickstart)
endif()
include(feelpp.module.pyfeelpp)
add_subdirectory(pyfeelpp)


###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 24 Jul 2018
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


if (FEELPP_ENABLE_PYTHON_WRAPPING )


set( SRC  pyexpr.cpp  )


file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelpython "${SRC};${HDR}")

install(PROGRAMS  ${CMAKE_CURRENT_SOURCE_DIR}/sympy2ginac.py   DESTINATION share/feelpp/data/python COMPONENT Bin)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/sympy2ginac.py DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp)
endif()

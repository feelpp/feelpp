###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Romain Hild <romain.hild@cemosis.fr>
#       Date: 09 Mar 2022
#
#  Copyright (C) 2022 Feel++ Consortium
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
cmake_minimum_required(VERSION 3.3)

message(STATUS "[pyfeelpp-mor] add pymodule pbdw")
feelpp_add_pymodule(pbdw SRCS pbdw.cpp DESTINATION feelpp/mor/pbdw LINK_LIBRARIES Feelpp::feelpp_mor )
add_dependencies(pyfeelppmor _pbdw)

install(FILES  __init__.py DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/feelpp/mor/pbdw)

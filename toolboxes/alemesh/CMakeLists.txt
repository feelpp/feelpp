###  CMakeLists.txt; coding: utf-8 --- 

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 18 Jan 2017
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


feelpp_add_application( alemesh_2d SRCS main.cpp
  LINK_LIBRARIES
  feelpp_modelcore
  DEFS FEELPP_DIM=2
  INSTALL
)

feelpp_add_application( alemesh_remesh_2d SRCS alemesh_remesh.cpp
  LINK_LIBRARIES
  feelpp_modelcore
  DEFS FEELPP_DIM=2
  INSTALL
  )

feelpp_add_application( alemesh_remesh_3d SRCS alemesh_remesh.cpp
  LINK_LIBRARIES
  feelpp_modelcore
  DEFS FEELPP_DIM=3
  INSTALL
  )


feelpp_add_testcase(cases PREFIX feelpp_toolbox_alemesh CATEGORY toolboxes/alemesh)


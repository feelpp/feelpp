###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
#       Date: 2013-07-22
#
#  Copyright (C) 2013
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

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
cmake_minimum_required(VERSION 2.8)
find_path(FEELPP_CMAKE_MODULES FindFeel++.cmake PATHS /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
if(FEELPP_CMAKE_MODULES)
	set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
else()
	message(FATAL_ERROR "Feel++ does not seem to have been installed on this platform")
endif()
find_package(Feel++)
endif()

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
feelpp_add_application(normalToBoundary SRCS normalToBoundary.cpp )
feelpp_add_application(dist2walls_simplex_2d SRCS dist2walls.cpp DEFS FEELPP_CONVEX=Simplex FEELPP_DIM=2 )
feelpp_add_application(dist2walls_hypercube_2d SRCS dist2walls.cpp DEFS FEELPP_CONVEX=Hypercube FEELPP_DIM=2 CFG dist2walls-hypercube.cfg)
feelpp_add_application(dist2wallsoptimized SRCS dist2wallsoptimized.cpp )
feelpp_add_application(dist2curve SRCS dist2curve.cpp CFG dist2curve.cfg )
feelpp_add_application(dist2curveperf SRCS dist2curveperf.cpp )
feelpp_add_application(dist2curveho SRCS dist2curveho.cpp )
feelpp_add_application(dist2curvehoperiodic SRCS dist2curvehoperiodic.cpp )
feelpp_add_application(dist2paramcurve SRCS dist2paramcurve.cpp )
feelpp_add_application(dist2paramsurface SRCS dist2paramsurface.cpp )

###  CMakeLists.txt; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2011-05-28
#
#  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
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
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/bamglib)

if(DRAWING)
  add_library(feel++_bamg
  Mesh2d.cpp
  bamg-gmsh.cpp
  bamglib/Mesh2.cpp
  bamglib/MeshGeom.cpp
  bamglib/MeshRead.cpp
  bamglib/Metric.cpp
  bamglib/R2.cpp
  bamglib/MeshQuad.cpp
  bamglib/MeshWrite.cpp
  bamglib/Meshio.cpp
  bamglib/QuadTree.cpp
  bamglib/SetOfE4.cpp
  bamglib/MeshDraw.cpp
  )
  else()
  add_library(feel++_bamg
  Mesh2d.cpp
  bamg-gmsh.cpp
  bamglib/Mesh2.cpp
  bamglib/MeshGeom.cpp
  bamglib/MeshRead.cpp
  bamglib/Metric.cpp
  bamglib/R2.cpp
  bamglib/MeshQuad.cpp
  bamglib/MeshWrite.cpp
  bamglib/Meshio.cpp
  bamglib/QuadTree.cpp
  bamglib/SetOfE4.cpp
  )
endif()


file(GLOB_RECURSE HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
INSTALL(FILES ${HDR}  DESTINATION include/feel/bamg COMPONENT Devel)



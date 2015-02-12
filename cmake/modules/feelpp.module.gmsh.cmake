###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2014-08-19
#
#  Copyright (C) 2014 Feel++ Consortium
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
if ( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gmsh )

  add_subdirectory(contrib/gmsh)
  # we include this directory : add some missing headers from Gmsh
  INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/contrib/gmsh )

#else( EXISTS ${CMAKE_SOURCE_DIR}/contrib/nlopt )
else( EXISTS ${CMAKE_SOURCE_DIR}/contrib/gmsh )
  FIND_PATH(GMSH_CONTRIB_INCLUDE_DIR BasisFactory.h
    $ENV{FEELPP_DIR}/include/feel/gmsh
    NO_DEFAULT_PATH)
  FIND_PATH(GMSH_CONTRIB_INCLUDE_DIR BasisFactory.h
    PATH_SUFFIXES feel/gmsh )
  include_directories(${GMSH_CONTRIB_INCLUDE_DIR})
  
endif(  EXISTS ${CMAKE_SOURCE_DIR}/contrib/gmsh  )
  

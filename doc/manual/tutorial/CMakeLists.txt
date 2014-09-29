###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2012-02-03
#
#  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
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
cmake_minimum_required (VERSION 2.8)


add_custom_target(
  loadmesh
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${FEELPP_SOURCE_DIR}/data/gmsh/feelpp/feel.geo
                                   ${FEELPP_BUILD_DIR}/doc/manual/tutorial
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(tutorial)

foreach(app myapp mymesh myintegrals myfunctionspace geotool loadmesh mylaplacian myadvection myexporter myexpression mymeshfromstring mybackend myfunctor )

  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg )
    feelpp_add_application( tut_${app} SRCS ${app}.cpp LABELS doc CFG ${app}.cfg)
  else()
    feelpp_add_application( tut_${app} SRCS ${app}.cpp LABELS doc)
  endif()
  install(TARGETS feelpp_tut_${app} DESTINATION bin/ COMPONENT Bin )

  add_dependencies(tutorial feelpp_tut_${app})
endforeach()

feelpp_add_application( tut_mystokes SRCS mystokes.cpp LABELS doc CFG mystokes.cfg bercovierengelmann.cfg kovasnay.cfg GEO mystokes.geo)
install(TARGETS feelpp_tut_mystokes DESTINATION bin/ COMPONENT Bin )
  add_dependencies(tutorial feelpp_tut_mystokes)


add_custom_target(tutorial-check COMMAND "ctest" ARGS "-VV" "-R" "tut")
add_dependencies(tutorial-check tutorial)
add_dependencies(check tutorial-check)

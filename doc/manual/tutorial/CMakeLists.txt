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

add_custom_target(
  loadmesh
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy ${FEELPP_SOURCE_DIR}/data/gmsh/feel.geo
                                   ${FEELPP_BUILD_DIR}/doc/manual/tutorial
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

foreach(app myapp mymesh myintegrals myfunctionspace geotool loadmesh mylaplacian myadvection myexporter myexpression )

  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg )
    feelpp_add_application( doc_${app} SRCS ${app}.cpp LABEL doc CFG ${app}.cfg)
  else()
    feelpp_add_application( doc_${app} SRCS ${app}.cpp LABEL doc)
  endif()
  install(TARGETS feelpp_doc_${app} DESTINATION bin/ COMPONENT Bin )
endforeach()
feelpp_add_application( doc_mystokes SRCS mystokes.cpp LABEL doc CFG mystokes.cfg bercovierengelmann.cfg kovasnay.cfg GEO mystokes.geo)
install(TARGETS feelpp_doc_mystokes DESTINATION bin/ COMPONENT Bin )

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
  COMMAND ${CMAKE_COMMAND} -E copy ${FEELPP_SOURCE_DIR}/data/medit/Cylref.mesh
                                   ${FEELPP_BUILD_DIR}/doc/manual/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_docex(myapp)
add_docex(mymesh)
add_docex(loadmesh)
add_docex(myintegrals)
add_docex(geotool)
add_docex(myfunctionspace)

if( NOT FEELPP_MINIMAL_CONFIGURATION )
  add_docex(laplacian)
  add_docex(dar)
  #add_docex(stokes)
  feelpp_add_application( doc_stokes SRCS stokes.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )
  feelpp_add_application( doc_stokes_lm SRCS stokes.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} DEFS FEELPP_USE_LM )
  feelpp_add_application( doc_stokes_poiseuille SRCS stokes_poiseuille.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )
  feelpp_add_application( doc_stokes_poiseuille_lm SRCS stokes_poiseuille.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} DEFS FEELPP_USE_LM )

  feelpp_add_application( doc_stokes3d SRCS stokes3d.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} )
  feelpp_add_application( doc_stokes3d_lm SRCS stokes3d.cpp LINK_LIBRARIES ${FEELPP_LIBRARIES} DEFS FEELPP_USE_LM )
endif()

# stokes_good stokes_locking stokes_spurious
# bratu nonlinearpow mortar multi_levels threefields nirb)

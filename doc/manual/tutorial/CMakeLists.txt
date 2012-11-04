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

foreach(app myapp mymesh myintegrals myfunctionspace geotool loadmesh )
  feelpp_add_application( doc_${app} SRCS ${app}.cpp LABEL doc)
endforeach()


if( NOT FEELPP_MINIMAL_CONFIGURATION )

  feelpp_add_application( doc_laplacian SRCS laplacian.cpp LABEL doc CFG laplacian.cfg SCRIPTS laplacian.py ADD_OT )

  feelpp_add_application( doc_laplacian_periodic SRCS periodic.cpp LABEL doc  )
  feelpp_add_application( doc_stokes_periodic SRCS stokes_periodic.cpp LABEL doc  )

  feelpp_add_application( doc_laplacian_full SRCS laplacian_full.cpp LABEL doc)
  feelpp_add_application( doc_laplacian_lagrange_multiplier SRCS laplacian_lagrange_multiplier.cpp LABEL doc)
  feelpp_add_application( doc_sound SRCS sound.cpp room.cpp LABEL doc )
  feelpp_add_application( doc_wave SRCS wave.cpp room.cpp LABEL doc )
  feelpp_add_application( doc_advection SRCS advection.cpp CFG advection.cfg LABEL doc )

  feelpp_add_application( doc_stokes SRCS stokes.cpp LABEL doc )
  feelpp_add_application( doc_stokes_lm SRCS stokes.cpp LABEL doc DEFS FEELPP_USE_LM )

  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet SRCS stokes_Dirichlet_Dirichlet.cpp LABEL doc )
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_lm SRCS stokes_Dirichlet_Dirichlet.cpp LABEL doc DEFS FEELPP_USE_LM )
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_cylinder SRCS stokes_Dirichlet_Dirichlet.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=2)
  feelpp_add_application( doc_stokes_Dirichlet_Dirichlet_hypercube SRCS stokes_Dirichlet_Dirichlet.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=1)

  feelpp_add_application( doc_stokes_Dirichlet_Neumann SRCS stokes_Dirichlet_Neumann.cpp LABEL doc )
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_lm SRCS stokes_Dirichlet_Neumann.cpp LABEL doc DEFS FEELPP_USE_LM )
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_cylinder SRCS stokes_Dirichlet_Neumann.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=2)
  feelpp_add_application( doc_stokes_Dirichlet_Neumann_hypercube SRCS stokes_Dirichlet_Neumann.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=1)

  feelpp_add_application( doc_stokes_Neumann_Neumann SRCS stokes_Neumann_Neumann.cpp LABEL doc )
  feelpp_add_application( doc_stokes_Neumann_Neumann_aorta SRCS stokes_Neumann_Neumann.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=3)
  feelpp_add_application( doc_stokes_Neumann_Neumann_cylinder SRCS stokes_Neumann_Neumann.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=2)
  feelpp_add_application( doc_stokes_Neumann_Neumann_hypercube SRCS stokes_Neumann_Neumann.cpp LABEL doc DEFS STOKESPRESSMESHTYPE=1)

  feelpp_add_application( doc_stokes_laplacian_inlet SRCS stokes_laplacian_inlet.cpp LABEL doc)
endif()


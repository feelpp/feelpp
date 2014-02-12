###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2013-02-18
#
#  Copyright (C) 2013 Feel++ Consortium
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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

feelpp_add_application( doc_laplacian SRCS laplacian.cpp LABEL doc CFG laplacian.cfg SCRIPTS laplacian.py ADD_OT )#
feelpp_add_application( doc_laplacian_full SRCS laplacian_full.cpp LABEL doc NO_MPI_TEST )
feelpp_add_application( doc_laplacian_mixed SRCS laplacian_mixed.cpp LABEL doc )
feelpp_add_application( doc_laplacian_lagrange_multiplier SRCS laplacian_lagrange_multiplier.cpp LABEL doc)
feelpp_add_application( doc_laplacian_periodic SRCS laplacian_periodic.cpp LABEL doc NO_MPI_TEST )
feelpp_add_application( doc_laplacian_harmonic SRCS harmonic.cpp LABEL doc  )

feelpp_add_application( doc_laplacian_quad SRCS laplacian_quad.cpp LABEL doc  CFG laplacian_polar.cfg)
feelpp_add_application( doc_laplacian_with_holes SRCS laplacian_with_holes.cpp LABEL doc  )

# DG examples
feelpp_add_application( doc_laplacian_dg SRCS laplacian_dg.cpp CFG laplacian_dg.cfg LABEL doc  )
feelpp_add_application( doc_laplacian_dg2 SRCS laplacian_dg2.cpp CFG laplacian_dg2.cfg LABEL doc  )

feelpp_add_application( doc_laplacian_polar SRCS laplacian_polar.cpp LABEL doc GEO polar.geo  CFG laplacian_polar.cfg DEFS FEELPP_POLAR)
feelpp_add_application( doc_laplacian_cartesian SRCS laplacian_polar.cpp LABEL doc GEO cartesian.geo CFG laplacian_cartesian.cfg )

feelpp_add_application( doc_laplacian_matlab SRCS laplacian_matlab.cpp GEO feel.geo SCRIPTS laplacian_matlab.m )

feelpp_add_application( doc_laplacian_block SRCS laplacian_block.cpp LABEL doc CFG laplacian_block.cfg )
feelpp_add_application( doc_laplacian_parabolic SRCS laplacian_parabolic.cpp LABEL doc CFG laplacian_parabolic.cfg GEO laplacian_parabolic.geo )

feelpp_add_application( doc_laplacian_substructured SRCS laplacian_substructured.cpp LABEL doc  )
feelpp_add_application( doc_laplacian_hypercube SRCS laplacian_hypercube.cpp LABEL doc CFG laplacian_hypercube.cfg )

feelpp_add_application( doc_laplacian_dirac SRCS laplacian_dirac.cpp LABEL doc CFG laplacian_dirac.cfg )

feelpp_add_application( doc_laplacian_lagrange_multiplier2 SRCS laplacian_lagrange_multiplier2.cpp LABEL doc CFG laplacian_lagrange_multiplier2.cfg)

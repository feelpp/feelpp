###  TEMPLATE.txt.tpl; coding: utf-8 ---
#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2013-02-18
#
#  Copyright (C) 2013-2016 Feel++ Consortium
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

feelpp_add_application( laplacian SRCS laplacian.cpp LABELS doc CFG laplacian.cfg SCRIPTS laplacian.py ADD_OT )#
feelpp_add_application( laplacian_full SRCS laplacian_full.cpp LABELS doc NO_MPI_TEST )
feelpp_add_application( laplacian_mixed SRCS laplacian_mixed.cpp LABELS doc CFG laplacian_mixed.cfg)
feelpp_add_application( laplacian_lagrange_multiplier SRCS laplacian_lagrange_multiplier.cpp LABELS doc CFG laplacian_lagrange_multiplier.cfg )
feelpp_add_application( laplacian_periodic SRCS laplacian_periodic.cpp LABELS doc NO_MPI_TEST )
feelpp_add_application( laplacian_harmonic SRCS harmonic.cpp LABELS doc  )

feelpp_add_application( laplacian_quad SRCS laplacian_quad.cpp LABELS doc  CFG laplacian_quad.cfg)
feelpp_add_application( laplacian_with_holes SRCS laplacian_with_holes.cpp LABELS doc  )

# DG examples
feelpp_add_application( laplacian_dg SRCS laplacian_dg.cpp CFG laplacian_dg.cfg LABELS doc  )
feelpp_add_application( laplacian_dg2 SRCS laplacian_dg2.cpp CFG laplacian_dg2.cfg LABELS doc  )

feelpp_add_application( laplacian_polar SRCS laplacian_polar.cpp LABELS doc GEO polar.geo  CFG laplacian_polar.cfg DEFS FEELPP_POLAR)
feelpp_add_application( laplacian_cartesian SRCS laplacian_polar.cpp LABELS doc GEO cartesian.geo CFG laplacian_cartesian.cfg )

feelpp_add_application( laplacian_matlab SRCS laplacian_matlab.cpp GEO feel.geo SCRIPTS laplacian_matlab.m )

feelpp_add_application( laplacian_block SRCS laplacian_block.cpp LABELS doc CFG laplacian_block.cfg )
feelpp_add_application( laplacian_parabolic SRCS laplacian_parabolic.cpp LABELS doc CFG laplacian_parabolic.cfg GEO laplacian_parabolic.geo )

feelpp_add_application( laplacian_hypercube SRCS laplacian_hypercube.cpp LABELS doc CFG laplacian_hypercube.cfg )

feelpp_add_application( laplacian_submesh SRCS laplacian_submesh.cpp LABELS doc  )

feelpp_add_application( laplacian_dirac SRCS laplacian_dirac.cpp LABELS doc CFG laplacian_dirac.cfg )

feelpp_add_application( laplacian_lagrange_multiplier2 SRCS laplacian_lagrange_multiplier2.cpp LABELS doc CFG laplacian_lagrange_multiplier2.cfg)

feelpp_add_application( laplacian_bean SRCS laplacian_bean.cpp LABELS doc CFG laplacian_bean.cfg)

feelpp_add_application( laplacian_meshsupport SRCS laplacian_meshsupport.cpp GEO laplacian_meshsupport.geo CFG laplacian_meshsupport.cfg )

# why this example here???
feelpp_add_application( laplacian_partition SRCS partition.cpp LABELS doc CFG partition.cfg)

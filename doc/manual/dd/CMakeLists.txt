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
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


feelpp_add_application( doc_traces SRCS traces.cpp LABELS doc   )
feelpp_add_application( doc_partitioning SRCS partitioning.cpp LABELS doc   )
feelpp_add_application( doc_explicit SRCS explicit.cpp LABELS doc   )


# Overlapping schwartz domain decomposition method
feelpp_add_application( doc_overlapping SRCS overlapping_ddmethodgeo.cpp overlapping_ddmethod.cpp LABELS doc NO_MPI_TEST  )

# Non-overlapping schwartz domain decomposition method
feelpp_add_application( doc_nonoverlapping SRCS nonoverlappinggeo.cpp nonoverlapping.cpp LABELS doc NO_MPI_TEST  )

# eigenmodes of Dirichlet to Neumann operator
feelpp_add_application( doc_multi_levels SRCS multi_levels.cpp LABELS doc NO_MPI_TEST  )

# Non-conforming domain decomposition methods
feelpp_add_application( doc_mortar SRCS mortar.cpp LABELS doc  NO_MPI_TEST )
feelpp_add_application( doc_threefields SRCS threefields.cpp TEST -pc_type ilu -ksp_monitor  LABELS doc NO_MPI_TEST )

if(FEELPP_HAS_HPDDM)
  feelpp_add_application(doc_geneo SRCS geneo.cpp geneo_coefficients.cpp CFG geneo.cfg LABELS doc)
endif()

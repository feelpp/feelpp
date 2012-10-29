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
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Overlapping schwartz domain decomposition method
feelpp_add_application( doc_overlapping SRCS overlapping_ddmethodgeo.cpp overlapping_ddmethod.cpp LABEL doc  )

# Non-overlapping schwartz domain decomposition method
feelpp_add_application( doc_nonoverlapping SRCS nonoverlappinggeo.cpp nonoverlapping.cpp LABEL doc  )

# eigenmodes of Dirichlet to Neumann operator
feelpp_add_application( doc_multi_levels SRCS multi_levels.cpp LABEL doc  )

# Non-conforming domain decomposition methods
feelpp_add_application( doc_mortar SRCS mortar.cpp LABEL doc  )
feelpp_add_application( doc_threefields SRCS threefields.cpp NO_TEST LABEL doc )

add_test( threefields-1 feelpp_doc_threefields -pc_type ilu -ksp_monitor )
set_property(TEST threefields-1 PROPERTY LABEL doc)

add_dependencies( doc
  feelpp_doc_overlapping
  feelpp_doc_nonoverlapping
  feelpp_doc_multi_levels
  feelpp_doc_mortar
  feelpp_doc_threefields )

###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2012-05-23
#
#  Copyright (C) 2012 Universit� Joseph Fourier (Grenoble I)
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
set_directory_properties(PROPERTIES LABEL testmor )

foreach( geo test_deim.geo test_geim.geo parameterspace.json )
  add_custom_target(copy-${geo})
  add_custom_command(
        TARGET copy-${geo} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/${geo}
                ${CMAKE_CURRENT_BINARY_DIR}/${geo})
endforeach()

feelpp_add_test( eim NO_MPI_TEST LINK_LIBRARIES Feelpp::feelpp_mor )
feelpp_add_test( eim_grepl LINK_LIBRARIES Feelpp::feelpp_mor CLI "--eim.dimension-max=15 --trainset-eim-size=10 --eim.rebuild-database=true" )
feelpp_add_test( eim_composite CFG test_eim_composite.cfg LINK_LIBRARIES Feelpp::feelpp_mor NO_MPI_TEST )

feelpp_add_test( deim   LINK_LIBRARIES Feelpp::feelpp_mor NO_MPI_TEST CLI "--gmsh.filename ${CMAKE_CURRENT_SOURCE_DIR}/test_deim.geo")

feelpp_add_test( rbspace LINK_LIBRARIES Feelpp::feelpp_mor  )
feelpp_add_test( rbspacevector LINK_LIBRARIES Feelpp::feelpp_mor )

feelpp_add_test( db CFG test_db.cfg LINK_LIBRARIES Feelpp::feelpp_mor CLI "-st_pc_factor_mat_solver_package mumps --gmsh.filename ${CMAKE_CURRENT_SOURCE_DIR}/test_db.geo" )

feelpp_add_test( parameterspace LINK_LIBRARIES Feelpp::feelpp_mor CLI "--config-file ${CMAKE_CURRENT_SOURCE_DIR}/parameterspace.cfg" )
feelpp_add_test( reducedglobaldof LINK_LIBRARIES Feelpp::feelpp_mor )

feelpp_add_test( geim LINK_LIBRARIES Feelpp::feelpp_mor CLI "--gmsh.filename ${CMAKE_CURRENT_SOURCE_DIR}/test_geim.geo")
feelpp_add_test( pbdw LINK_LIBRARIES Feelpp::feelpp_mor CLI "--gmsh.filename ${CMAKE_CURRENT_SOURCE_DIR}/test_geim.geo")

if (FEELPP_HAS_GLPK)
  feelpp_add_test( empiricalquadrature LINK_LIBRARIES Feelpp::feelpp_mor )
endif()
feelpp_add_test( exprevaluator LINK_LIBRARIES Feelpp::feelpp_mor )

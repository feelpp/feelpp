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
set_directory_properties(PROPERTIES LABEL testcrb )
feelpp_add_test( eim NO_MPI_TEST )
feelpp_add_test( eim_grepl CLI "--eim.dimension-max=15 --trainset-eim-size=10 --eim.rebuild-database=true" )
feelpp_add_test( eim_composite CFG test_eim_composite.cfg NO_MPI_TEST )

feelpp_add_test( deim GEO test_deim.geo NO_MPI_TEST)

feelpp_add_test( rbspace )
feelpp_add_test( rbspacevector )

feelpp_add_test( db CFG test_db.cfg GEO test_db.geo CLI "-st_pc_factor_mat_solver_package mumps" )

feelpp_add_test( parameterspace )
feelpp_add_test( reducedglobaldof )

feelpp_add_test( geim GEO test_geim.geo )
feelpp_add_test( pbdw GEO test_geim.geo )

if (FEELPP_HAS_GLPK)
  feelpp_add_test( empiricalquadrature )
endif()
feelpp_add_test( exprevaluator )
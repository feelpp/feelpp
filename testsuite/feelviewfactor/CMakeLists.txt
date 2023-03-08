# ##  TEMPLATE.txt.tpl; coding: utf-8 ---

# Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# Date: 2022-06-08
#
# Copyright (C) 2022 Feel++ Consortium
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
set_directory_properties(PROPERTIES LABEL testviewfactor)
feelpp_add_test(viewfactor_quadrature NO_MPI_TEST CFG cases/viewfactor.cfg )
feelpp_add_testcase(cases PREFIX testviewfactors)
add_dependencies(feelpp_test_viewfactor_quadrature testviewfactors_add_testcase_cases )

set_directory_properties(PROPERTIES LABEL testviewfactor)
feelpp_add_test(viewfactor_raytracing NO_MPI_TEST CFG cases/viewfactor.cfg )
add_dependencies(feelpp_test_viewfactor_raytracing testviewfactors_add_testcase_cases )
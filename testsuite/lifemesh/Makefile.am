# -*- mode: makefile -*-
#
#  This file is part of the LifeV library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#       Date: 2005-02-01
#
#  Copyright (C) 2005 EPFL
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
include $(top_srcdir)/admin/common.am

AM_LDFLAGS=-static
check_PROGRAMS  = test_entity test_mesh test_regiontree test_mesh_codim1

test_entity_SOURCES = test_entity.cpp
test_entity_LDADD   = ../../life/lifemesh/liblifemesh.la ../../life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_mesh_SOURCES = test_mesh.cpp
test_mesh_LDADD   =  ../../life/liblife.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_mesh_codim1_SOURCES = test_mesh_codim1.cpp
test_mesh_codim1_LDADD   =  ../../life/liblife.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_regiontree_SOURCES = test_regiontree.cpp
test_regiontree_LDADD   = $(life_libs)



EXTRA_DIST = testsuite.at

include $(top_srcdir)/testsuite/Makefile.testsuite

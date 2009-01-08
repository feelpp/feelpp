# -*- makefile -*-
#
# Author: Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#
# Copyright (C) 2004 EPFL
#
#  Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
include $(top_srcdir)/admin/common.am
include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc

check_PROGRAMS 	= test_mpi test_boostmpi test_traits test_singleton test_factory test_debug test_openmp test_context test_material

test_mpi_SOURCES    = test_mpi.cpp
test_mpi_LDADD   = $(top_builddir)/life/lifecore/liblifecore.la  $(petsc_libs)
test_boostmpi_SOURCES    = test_boostmpi.cpp
test_boostmpi_LDADD   = $(top_builddir)/life/lifecore/liblifecore.la  $(petsc_libs)

test_traits_SOURCES = test_traits.cpp
test_traits_LDADD   = $(top_builddir)/life/lifecore/liblifecore.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)


test_debug_SOURCES  = test_debug.cpp

test_singleton_SOURCES	= test_singleton.cpp
test_singleton_LDADD   = $(top_builddir)/life/lifecore/liblifecore.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_factory_SOURCES	= test_factory.cpp
test_openmp_SOURCES	= test_openmp.cpp

test_context_SOURCES	= test_context.cpp
test_context_LDADD   = $(top_builddir)/life/lifecore/liblifecore.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)


test_material_SOURCES	= test_material.cpp
test_material_LDADD   = $(top_builddir)/life/lifematerial/liblifematerial.la  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)



EXTRA_DIST			= testsuite.at



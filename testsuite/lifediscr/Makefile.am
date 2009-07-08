# -*- makefile -*-
#
#
#
# Copyright (C) 2004 Christophe Prud'homme (christophe.prudhomme@epfl.ch)
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


# AM_CPPFLAGS         = -I../../life/lifefem -I$(top_srcdir)/boost
noinst_HEADERS		=

EXTRA_DIST =  testsuite.at test_integration_ho_geometry.geo

check_HEADERS       =
check_PROGRAMS		= test_bdf \
                      test_interpolation \
                      test_integration test_integration_ho test_integration_3d \
                      test_space \
                      test_geomap test_mixed


test_mixed_SOURCES	= test_mixed.cpp
test_mixed_LDADD      = ../../life/liblife.la

test_bdf_SOURCES	= test_bdf.cpp
test_bdf_LDADD      = ../../life/liblife.la
test_bdf_LDFLAGS    = -static

test_space_SOURCES	= test_space.cpp
test_space_LDADD    = $(life_libs) $(boost_unit_test_framework_lib)  $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_integration_SOURCES	= test_integration.cpp
test_integration_LDADD    = ../../life/liblife.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

#test_integration_dirac_SOURCES	= test_integration_dirac.cpp
#test_integration_dirac_LDADD    = ../../life/liblife.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_integration_ho_SOURCES	= test_integration_ho.cpp
test_integration_ho_LDADD    = ../../life/liblife.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_interpolation_SOURCES	= test_interpolation.cpp
test_interpolation_LDADD    = ../../life/liblife.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_integration_3d_SOURCES	= test_integration_3d.cpp
test_integration_3d_LDADD    = ../../life/liblife.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_geomap_SOURCES     = test_geomap.cpp
test_geomap_LDADD       = $(life_libs)


# test_fem_SOURCES	= test_fem.cpp
# test_fem_LDADD		= $(life_libs) $(boost_libs) $(boost_unit_test_framework_lib)


# test_spectral_convergence_2D_SOURCE = test_spectral_convergence_2D.cpp
# test_spectral_convergence_2D_LDADD = $(life_libs) $(boost_libs) $(boost_unit_test_framework_lib)

# test_spectral_convergence_3D_SOURCE = test_spectral_convergence_3D.cpp
# test_spectral_convergence_3D_LDADD = $(life_libs) $(boost_libs) $(boost_unit_test_framework_lib)

# test_lag_precond_SOURCE = test_lag_precond.cpp
# test_lag_precond_LDADD = $(life_libs) $(boost_libs) $(boost_unit_test_framework_lib)

# test_mult_petsc_SOURCE = test_mult_petsc.cpp
# test_mult_petsc_LDADD = $(life_libs) $(parmetis_libs) $(boost_libs) #$(boost_unit_test_framework_lib)

# test_prosi_SOURCE = test_prosi.cpp
# test_prosi_LDADD = $(life_libs) $(boost_libs)

# #AM_CXXFLAGS = -DBOOST_UBLAS_NO_ELEMENT_PROXIES -DNDEBUG

# POTFILES	= #data



check-local:
	@if test "$(top_srcdir)" != "$(top_builddir)"; then \
     for i in $(srcdir)/*.geo $(srcdir)/*.py $(srcdir)/*.m; do \
	 ln -sf $$i;\
     done;\
     fi


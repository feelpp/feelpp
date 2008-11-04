# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#       Date: 2005-07-28
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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc

AM_CPPFLAGS         = -I../../life/lifepoly -I$(top_srcdir)/boost #-I/usr/include/vtk/

check_HEADERS       = testlagrange.hpp

check_PROGRAMS		= test_jacobi test_imsimplex test_im test_poly test_qm_qk \
                             test_qm_pk test_qm_qk_3D test_qm_pk_3D   \
			  test_qm_pk_double test_lag
#			  test_spectral_1D test_spectral_2D test_spectral_3D test_spectral_laplacian_3D

test_jacobi_SOURCES	= test_jacobi.cpp
test_jacobi_LDADD   =   $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la   $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_im_SOURCES	= test_im.cpp
test_im_LDADD   = $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(boost_filesystem_lib) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_imsimplex_SOURCES	= test_imsimplex.cpp
test_imsimplex_LDADD   = $(top_builddir)/life/liblife.la
test_imsimplex_LDFLAGS=-static

test_poly_SOURCES       = test_poly.cpp
test_poly_LDADD   = $(top_builddir)/life/lifefilters/liblifefilters.la $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(boost_filesystem_lib) $(boost_serialization_lib) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) #$(boost_unit_test_framework_lib)

include $(top_srcdir)/testsuite/Makefile.testsuite

AM_CXXFLAGS = -DBOOST_UBLAS_NO_ELEMENT_PROXIES #-DNDEBUG

test_qm_qk_SOURCES	= test_qm_qk.cpp
test_qm_qk_LDADD   = $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_qm_pk_SOURCES	= test_qm_pk.cpp
test_qm_pk_LDADD   = $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_qm_qk_3D_SOURCES      = test_qm_qk_3D.cpp
test_qm_qk_3D_LDADD   =  $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_qm_pk_3D_SOURCES      = test_qm_pk_3D.cpp
test_qm_pk_3D_LDADD   = $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_qm_pk_double_SOURCES	= test_qm_pk_double.cpp
test_qm_pk_double_LDADD   =  $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_lag_SOURCES = test_lag.cpp
test_lag_LDADD   =   $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

# test_spectral_1D_SOURCES = test_spectral_1D.cpp
# test_spectral_1D_LDADD   =  $(life_libs) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

# test_spectral_2D_SOURCES = test_spectral_2D.cpp
# test_spectral_2D_LDADD   =  $(life_libs) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

# test_spectral_3D_SOURCES = test_spectral_3D.cpp
# test_spectral_3D_LDADD   =  $(life_libs) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

# test_spectral_laplacian_3D_SOURCES = test_spectral_laplacian_3D.cpp
# test_spectral_laplacian_3D_LDADD   =  $(life_libs) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

#test_dubiner_SOURCES    = test_dubiner.cpp
#test_dubiner_LDADD   =   $(top_builddir)/life/lifemesh/liblifemesh.la $(top_builddir)/life/lifecore/liblifecore.la $(boost_filesystem_lib) #$(boost_unit_test_framework_lib)

EXTRA_DIST			= testsuite.at
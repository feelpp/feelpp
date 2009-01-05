# -*- makefile -*-
#
# Author: Christophe Prud'homme
#
# Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano
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
include $(top_srcdir)/admin/common.am

LIBS=

check_HEADERS		= MatrixTest.hpp

check_LTLIBRARIES	= liblifealg_testsuite.la

LDADD= $(top_builddir)/life/liblife.la
liblifealg_testsuite_la_SOURCES	=
#liblifealg_testsuite_la_LIBADD	=

check_PROGRAMS	= test_petsc  test_umfpack  test_ublas test_matrix

test_petsc_SOURCES	= test_petsc.cpp
test_petsc_LDADD    = $(LDADD) $(boost_libs) $(boost_unit_test_framework_lib)

test_umfpack_SOURCES	= test_umfpack.cpp
test_umfpack_LDADD    = $(LDADD) $(boost_libs) $(boost_unit_test_framework_lib)

test_ublas_SOURCES          = test_ublas.cpp
test_ublas_LDADD    = $(LDADD) $(top_builddir)/boost/libs/log/src/libboost_log.la -lpthread -llapack -lblas $(boost_filesystem_libs)

test_matrix_SOURCES=test_matrix.cpp
test_matrix_LDADD=$(LDADD) $(boost_filesystem_lib) $(petsc_libs)

#-lmkl_solver -lmkl_lapack -lmkl -lmkl_p4 -lmkl_ia32 -lguide /usr/lib/libg2c.so




EXTRA_DIST  = testsuite.at

MOSTLYCLEANFILES=mass.m convdiff.m

# include $(top_srcdir)/testsuite/Makefile.testsuite

# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-05-08
#
#  Copyright (C) 2009 Université Joseph Fourier
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
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
include $(top_srcdir)/admin/common.am
include $(top_srcdir)/testsuite/Makefile.testsuite

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc

check_PROGRAMS 	= test_material

test_material_SOURCES    = test_material.cpp
test_material_LDADD   = $(LIFE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

EXTRA_DIST			= testsuite.at




# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2006-12-30
#
#  Copyright (C) 2006 Université Joseph Fourier
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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc


check_PROGRAMS 	= test_importergmsh

test_importergmsh_SOURCES= test_importergmsh.cpp
test_importergmsh_LDADD=$(life_libs)

EXTRA_DIST = testsuite.at

MOSTLYCLEANFILES=

include $(top_srcdir)/testsuite/Makefile.testsuite

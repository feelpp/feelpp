# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2008-11-11
#
#  Copyright (C) 2008 Université Joseph Fourier
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

ACLOCAL_AMFLAGS = -I admin

bin_PROGRAMS=myapp mymesh myintegrals laplacian

laplacian_SOURCES=laplacian.cpp
laplacian_LDADD=$(LIFE_LIBS)

myintegrals_SOURCES=myintegrals.cpp
myintegrals_LDADD=$(LIFE_LIBS)

mymesh_SOURCES=mymesh.cpp
mymesh_LDADD=$(LIFE_LIBS)

myapp_SOURCES=myapp.cpp
myapp_LDADD=$(LIFE_LIBS)


# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2008-05-25
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

lib_LTLIBRARIES		= liblifematerial.la

liblifematerial_la_SOURCES	=  material.cpp materiallib.cpp
liblifematerial_la_LDFLAGS  = -no-undefined
liblifematerial_la_LIBADD   = ../lifecore/liblifecore.la -lboost_serialization

pkglib_LTLIBRARIES = stdmaterial.la


stdmaterial_la_SOURCES = air.cpp castiron.cpp
stdmaterial_la_LDFLAGS = -module
stdmaterial_la_LIBADD = liblifematerial.la -lboost_serialization


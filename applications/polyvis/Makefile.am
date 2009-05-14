# -*- mode: makefile -*-
#
#  This file is part of the LifeV library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#       Date: 2006-02-01
#
#  Copyright (C) 2006 EPFL
#  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)
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

ACLOCAL_AMFLAGS = -I admin

bin_PROGRAMS    = life_polyvis

life_polyvis_LDADD    = $(LIFE_LIBS)
life_polyvis_SOURCES  = polyvisapp.cpp polyvisbase.cpp \
                        lagrange1d.cpp lagrange2d.cpp lagrange3d.cpp\
                        orthonormalpset1d.cpp orthonormalpset2d.cpp orthonormalpset3d.cpp


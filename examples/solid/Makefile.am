# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-02-18
#
#  Copyright (C) 2007 Université Joseph Fourier
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

check_PROGRAMS=beam stvenant_kirchhoff axi1D # stvenant_kirchhoff_ho

#CXXFLAGS= -DLIFE_USE_BOOST_LAMBDA

beam_SOURCES=beam.cpp
beam_LDADD = $(LIFE_LIBS) -lboost_signals

axi1D_SOURCES=axi_1D.cpp
axi1D_LDADD = $(LIFE_LIBS) -lboost_signals

stvenant_kirchhoff_SOURCES=stvenant_kirchhoff.cpp
stvenant_kirchhoff_LDADD = $(LIFE_LIBS) -lboost_signals

#stvenant_kirchhoff_ho_SOURCES=stvenant_kirchhoff_ho.cpp
#stvenant_kirchhoff_ho_LDADD = $(LIFE_LIBS) -lboost_signals

EXTRA_DIST=solid.at

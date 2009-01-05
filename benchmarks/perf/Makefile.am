# -*- mode: makefile -*-
#
#  This file is part of the LifeV library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-11-16
#
#  Copyright (C) 2007 Université Joseph Fourier Grenoble 1
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
#AM_CPPFLAGS =  $(AM_CPPFLAGS) -I$(srcdir)/../

EXTRA_DIST=perf.at


check_PROGRAMS    = perf
perf_SOURCES  = bench1.cpp
perf_LDADD    = ../liblifelog.la $(LIFE_LIBS)  -lpthread
perf_LDFLAGS  = -static

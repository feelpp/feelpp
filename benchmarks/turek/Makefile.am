# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2008-05-01
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

check_PROGRAMS=turek

#AM_CXXFLAGS=-ftree-vectorizer-verbose=1 -ftree-vectorize
turek_SOURCES=turekapp.cpp data.cpp \
              turek_d2o2g1.cpp

#turek_d2o2g2.cpp \
#              turek_d2o3g1.cpp turek_d2o3g2.cpp \
3              turek_d2o4g1.cpp turek_d2o4g2.cpp \
3              turek_d3o2g1.cpp turek_d3o2g2.cpp \
3              turek_d3o3g1.cpp turek_d3o3g2.cpp \
3              turek_d3o4g1.cpp turek_d3o4g2.cpp

turek_LDADD = ../../life/liblife.la -lprofiler
turek_LDFLAGS = -static -lprofiler



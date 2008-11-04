# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-02-23
#
#  Copyright (C) 2007 Universit� Joseph Fourier
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

SUBDIRS=. convergence perf kovasznay stokes turek structure

lib_LTLIBRARIES=liblifelog.la

check_HEADERS=logs.hpp
liblifelog_la_SOURCES=logs.cpp
liblifelog_la_LIBADD=../boost/libs/log/src/libboost_log.la

EXTRA_DIST=benchmarksuite.at
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
include $(top_srcdir)/admin/common.am
check_PROGRAMS=convection splitting2d cavity2d rclux #test_nsale test_ale1  test_ale2  test_ale3  test_ale4 # twofluid2d cavity2d cavity3d

noinst_HEADERS=advreact.hpp oseen.hpp  splitting.hpp  twofluid.hpp

convection_SOURCES=convection.cpp
convection_LDADD = $(LIFE_LIBS)

cavity2d_SOURCES=cavity2d.cpp
cavity2d_LDADD = $(LIFE_LIBS)

#cavity3d_SOURCES=cavity3d.cpp
#cavity3d_LDADD = $(life_libs)

splitting2d_SOURCES=splitting2d.cpp
splitting2d_LDADD = $(LIFE_LIBS)

rclux_SOURCES=rcluxmesh.cpp rclux.cpp
rclux_LDADD = $(LIFE_LIBS)

#test_nsale_SOURCES=test_nsALE.cpp
#test_nsale_LDADD = $(LIFE_LIBS)


# test_ale1_SOURCES=test_ale1.cpp
# test_ale1_LDADD = $(LIFE_LIBS)

# test_ale2_SOURCES=test_ale2.cpp
# test_ale2_LDADD = $(LIFE_LIBS)


# test_ale3_SOURCES=test_ale3.cpp
# test_ale3_LDADD = $(LIFE_LIBS)


# test_ale4_SOURCES=test_ale4.cpp
# test_ale4_LDADD = $(LIFE_LIBS)

#twofluid2d_SOURCES=twofluid2d.cpp
#twofluid2d_LDADD = $(life_libs)

EXTRA_DIST=fluid.at
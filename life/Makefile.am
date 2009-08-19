# -*- Mode : makefile -*-
#
#    AUTHOR: Christophe Prud'homme
#       ORG: EPFL
#    E-MAIL: christophe.prudhomme@epfl.ch
#
# DESCRIPTION:
# ============
#  Distributed under the GPL(GNU Public License):
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# DESCRIP-END.
#
include $(top_srcdir)/admin/common.am

SUBDIRS		= lifecore  lifealg lifemesh lifepoly lifediscr lifevf  lifefilters lifesystem lifematerial


dist_pkgdata_DATA  = debug.areas

lib_LTLIBRARIES		= liblife.la

liblife_la_SOURCES	=  options.cpp
liblife_la_LIBADD   = lifematerial/liblifematerial.la \
                      lifefilters/liblifefilters.la \
                      lifediscr/liblifediscr.la     \
                      lifemesh/liblifemesh.la       \
                      lifealg/liblifealg.la         \
                      lifecore/liblifecore.la
liblife_la_LDFLAGS  = -no-undefined -version-info 0:0;0

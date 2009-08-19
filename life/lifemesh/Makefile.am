# -*- Mode : makefile -*-
#
# Copyright (C) 2004 EPFL, INRIA and Politechnico di Milano
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

include $(top_srcdir)/admin/common.am

#INCLUDES=-I/usr/include/metis
lib_LTLIBRARIES		= liblifemesh.la
liblifemesh_la_SOURCES	= marker.cpp bareitems.cpp \
                          entities.cpp simplex.cpp simplexproduct.cpp \
                          regiontree.cpp kdtree.cpp\
                          functors.cpp \
                          meshbase.cpp

liblifemesh_la_LIBADD	= ../lifecore/liblifecore.la
liblifemesh_la_LDFLAGS  = -version-info 0:0:0


#bin_PROGRAMS    = tt
#tt_SOURCES = tt.cpp
#tt_LDADD	= ../lifemesh/liblifemesh.la ../lifecore/liblifecore.la

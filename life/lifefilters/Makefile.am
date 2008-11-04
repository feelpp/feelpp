# -*- makefile -*-
#
# Author: Christophe Prud'homme
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

SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc


lib_LTLIBRARIES	= liblifefilters.la

liblifefilters_la_SOURCES	= gmsh.cpp gmshtensorizeddomain.cpp gmshsimplexdomain.cpp\
                              importergambit.cpp\
                              exporter.cpp exporterensight.cpp exportergmsh.cpp exportergnuplot.cpp


liblifefilters_la_LIBADD	=  $(BOOST_FILESYSTEM_LIB) ../../boost/libs/serialization/src/libboost_serialization.la ../../boost/libs/mpi/src/libboost_mpi.la


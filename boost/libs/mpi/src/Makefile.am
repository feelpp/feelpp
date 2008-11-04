# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2007-06-28
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
lib_LTLIBRARIES=libboost_mpi.la

libboost_mpi_la_SOURCES=broadcast.cpp communicator.cpp computation_tree.cpp content_oarchive.cpp environment.cpp exception.cpp graph_communicator.cpp group.cpp intercommunicator.cpp mpi_datatype_cache.cpp mpi_datatype_oarchive.cpp packed_iarchive.cpp packed_oarchive.cpp packed_skeleton_iarchive.cpp packed_skeleton_oarchive.cpp point_to_point.cpp request.cpp text_skeleton_oarchive.cpp timer.cpp #text_skeleton_iarchive.cpp

LIBS=
libboost_mpi_la_LIBADD=../../serialization/src/libboost_serialization.la $(MPILIBS)

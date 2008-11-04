# -*- mode: makefile -*-
#
#  This file is part of the Life library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2008-03-27
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
AM_LDFLAGS=-static
check_PROGRAMS=laplacian1dP1 laplacian1dP5\
               laplacian1d2dP1 laplacian1d2dP5\
               laplacian2dP1 laplacian2dP5\
               laplacian2d3dP1 \
               laplacian2dP1h laplacian2dP5h\
               laplacian3dP1 laplacian3dP2 laplacian3dP3 laplacian3dP5\
               laplacian3dP1h\
               \
               laplacian2dP1v laplacian2dP4v

check_HEADERS=laplacian.hpp
laplacian1dP1_SOURCES=laplacian1dP1.cpp
laplacian1dP1_LDADD=../../life/liblife.la

laplacian1dP5_SOURCES=laplacian1dP5.cpp
laplacian1dP5_LDADD=../../life/liblife.la

laplacian1d2dP1_SOURCES=laplacian1d2dP1.cpp
laplacian1d2dP1_LDADD=../../life/liblife.la

laplacian1d2dP5_SOURCES=laplacian1d2dP5.cpp
laplacian1d2dP5_LDADD=../../life/liblife.la

laplacian2dP1_SOURCES=laplacian2dP1.cpp
laplacian2dP1_LDADD=../../life/liblife.la

laplacian2dP5_SOURCES=laplacian2dP5.cpp
laplacian2dP5_LDADD=../../life/liblife.la

laplacian2d3dP1_SOURCES=laplacian2d3dP1.cpp
laplacian2d3dP1_LDADD=../../life/liblife.la

laplacian2dP1_SOURCES=laplacian2dP1.cpp
laplacian2dP1_LDADD=../../life/liblife.la

laplacian2dP1h_SOURCES=laplacian2dP1h.cpp
laplacian2dP1h_LDADD=../../life/liblife.la

laplacian2dP5h_SOURCES=laplacian2dP5h.cpp
laplacian2dP5h_LDADD=../../life/liblife.la

laplacian3dP1_SOURCES=laplacian3dP1.cpp
laplacian3dP1_LDADD=../../life/liblife.la

laplacian3dP2_SOURCES=laplacian3dP2.cpp
laplacian3dP2_LDADD=../../life/liblife.la

laplacian3dP3_SOURCES=laplacian3dP3.cpp
laplacian3dP3_LDADD=../../life/liblife.la

laplacian3dP5_SOURCES=laplacian3dP5.cpp
laplacian3dP5_LDADD=../../life/liblife.la

laplacian3dP1h_SOURCES=laplacian3dP1h.cpp
laplacian3dP1h_LDADD=../../life/liblife.la

laplacian2dP1v_SOURCES=laplacian2dP1v.cpp
laplacian2dP1v_LDADD=../../life/liblife.la

laplacian2dP4v_SOURCES=laplacian2dP4v.cpp
laplacian2dP4v_LDADD=../../life/liblife.la

EXTRA_DIST=convergence.at convergence.py
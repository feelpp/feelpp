# -*- Mode : makefile; indent-tabs-mode: t -*-
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

check_PROGRAMS=pmesh laplacian laplacian_p laplacian_p2 advection dar stokes wave laplacian_multiplicateur_lagrange periodic

# stokesT  test_ale  #test_formulation nonlinear
#CXXFLAGS= -DLIFE_USE_BOOST_LAMBDA




pmesh_SOURCES=pmesh.cpp
pmesh_LDADD =$(LIFE_LIBS)

laplacian_p_SOURCES=laplacian_p.cpp
laplacian_p_LDADD = $(LIFE_LIBS)

laplacian_p2_SOURCES=laplacian_p2.cpp
laplacian_p2_LDADD = $(LIFE_LIBS)

laplacian_SOURCES=laplacian.cpp
laplacian_LDADD = $(LIFE_LIBS)

laplacian_multiplicateur_lagrange_SOURCES=laplacian_multiplicateur_lagrange.cpp
laplacian_multiplicateur_lagrange_LDADD = $(LIFE_LIBS)
laplacian_multiplicateur_lagrange_LDFLAGS = -static

periodic_SOURCES=periodic.cpp
periodic_LDADD = $(LIFE_LIBS)
periodic_LDFLAGS = -static

advection_SOURCES=advection.cpp mesh.cpp
advection_LDADD = $(LIFE_LIBS)



dar_SOURCES=dar.cpp
dar_LDADD = $(LIFE_LIBS)

stokes_SOURCES=stokes.cpp
stokes_LDADD = $(LIFE_LIBS)
stokes_LDFLAGS=-static

#stokesT_SOURCES=stokesT.cpp
#stokesT_LDADD = $(LIFE_LIBS)

wave_SOURCES=wave.cpp
wave_LDADD = $(LIFE_LIBS)

if WITH_PETSC

check_PROGRAMS += bratu

bratu_SOURCES=bratu.cpp
bratu_LDADD = $(LIFE_LIBS)

endif

#test_ale_SOURCES=test_ale.cpp
#test_ale_LDADD = $(LIFE_LIBS)

#test_formulation_SOURCES=test_formulation.cpp
#test_formulation_LDADD = $(LIFE_LIBS)

#nonlinear_SOURCES=nonlinear.cpp
#nonlinear_LDADD = $(LIFE_LIBS)


if WITH_SLEPC

check_PROGRAMS += sound

sound_SOURCES=sound.cpp room.cpp
sound_LDADD = $(LIFE_LIBS)

endif



EXTRA_DIST=generic.at
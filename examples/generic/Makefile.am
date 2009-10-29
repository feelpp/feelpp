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

check_PROGRAMS=life_pmesh \
               life_laplacian \
			   life_laplacian_p\
			   life_laplacian_p2\
               life_advection \
               life_dar\
               life_stokes \
               life_wave \
               life_laplacian_multiplicateur_lagrange\
               life_periodic

# stokesT  test_ale  #test_formulation nonlinear
#CXXFLAGS= -DLIFE_USE_BOOST_LAMBDA




life_pmesh_SOURCES=pmesh.cpp
life_pmesh_LDADD =$(LIFE_LIBS)

life_laplacian_p_SOURCES=laplacian_p.cpp
life_laplacian_p_LDADD = $(LIFE_LIBS)

life_laplacian_p2_SOURCES=laplacian_p2.cpp
life_laplacian_p2_LDADD = $(LIFE_LIBS)

life_laplacian_SOURCES=laplacian.cpp
life_laplacian_LDADD = $(LIFE_LIBS)

life_laplacian_multiplicateur_lagrange_SOURCES=laplacian_multiplicateur_lagrange.cpp
life_laplacian_multiplicateur_lagrange_LDADD = $(LIFE_LIBS)
life_laplacian_multiplicateur_lagrange_LDFLAGS = -static

life_periodic_SOURCES=periodic.cpp
life_periodic_LDADD = $(LIFE_LIBS)
life_periodic_LDFLAGS = -static

life_advection_SOURCES=advection.cpp mesh.cpp
life_advection_LDADD = $(LIFE_LIBS)



life_dar_SOURCES=dar.cpp
life_dar_LDADD = $(LIFE_LIBS)

life_stokes_SOURCES=stokes.cpp
life_stokes_LDADD = $(LIFE_LIBS)
life_stokes_LDFLAGS=-static

#stokesT_SOURCES=stokesT.cpp
#stokesT_LDADD = $(LIFE_LIBS)

life_wave_SOURCES=wave.cpp
life_wave_LDADD = $(LIFE_LIBS)

if WITH_PETSC

check_PROGRAMS += life_bratu

life_bratu_SOURCES=bratu.cpp
life_bratu_LDADD = $(LIFE_LIBS)

endif

#test_ale_SOURCES=test_ale.cpp
#test_ale_LDADD = $(LIFE_LIBS)

#test_formulation_SOURCES=test_formulation.cpp
#test_formulation_LDADD = $(LIFE_LIBS)

#nonlinear_SOURCES=nonlinear.cpp
#nonlinear_LDADD = $(LIFE_LIBS)


if WITH_SLEPC

check_PROGRAMS += life_sound

life_sound_SOURCES=sound.cpp room.cpp
life_sound_LDADD = $(LIFE_LIBS)

endif



EXTRA_DIST=generic.at
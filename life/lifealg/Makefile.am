# -*- makefile -*-
#
# Author: Christophe Prud'homme
#
# Copyright (C) 2004 EPFL
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

if HAVE_PETSC
PETSC_SRCS=vectorpetsc.cpp matrixpetsc.cpp solverlinearpetsc.cpp solvernonlinearpetsc.cpp
endif
if HAVE_SLEPC
SLEPC_SRCS=solvereigenslepc.cpp
endif
if HAVE_TRILINOS
TRILINOS_SRCS=vectorepetra.cpp matrixepetra.cpp preconditionerml.cpp preconditionerifpack.cpp backendtrilinos.cpp solvernonlineartrilinos.cpp
endif
if HAVE_UMFPACK
UMFPACK_SRCS=solverumfpack.cpp
endif
#if HAVE_PARDISO
PARDISO_SRCS= #solverpardiso.cpp
#endif

lib_LTLIBRARIES	= liblifealg.la

liblifealg_la_SOURCES	= dummy.cpp datamap.cpp \
                          vector.cpp vectorublas.cpp \
                          backend.cpp backendgmm.cpp backend_adaptive_reuse_pc.cpp backendpetsc.cpp \
                          graphcsr.cpp matrixgmm.cpp solvereigen.cpp \
                          solvernonlinear.cpp \
                          $(UMFPACK_SRCS)\
                          $(PETSC_SRCS)\
                          $(SLEPC_SRCS)\
                          $(TRILINOS_SRCS)\
                          $(PARDISO_SRCS)

liblifealg_la_LIBADD	=  ../lifecore/liblifecore.la $(petsc_libs) $(umfpack_libs) $(aztec_libs)
liblifealg_la_LDFLAGS  = -version-info 0:0:0


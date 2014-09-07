# -*- mode: cmake; coding: utf-8 -*-
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2009-11-29
#
#  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
#
if( PETSC_FOUND )
  set( FEELPP_PETSC_SRCS vectorpetsc.cpp matrixpetsc.cpp solverlinearpetsc.cpp solvernonlinearpetsc.cpp functionspetsc.cpp preconditionerpetsc.cpp)

  if ( FEELPP_HAS_SLEPC )
    set( FEELPP_SLEPC_SRCS solvereigenslepc.cpp )
  endif()
endif()

if ( TRILINOS_LIBRARIES )
  set( FEELPP_TRILINOS_SRCS vectorepetra.cpp matrixepetra.cpp preconditionerml.cpp preconditionerifpack.cpp backendtrilinos.cpp solvernonlineartrilinos.cpp )
endif()

set( SRC
  dummy.cpp enums.cpp datamap.cpp
  vector.cpp vectorublas.cpp
  preconditioner.cpp
  backend.cpp  matrixshellsparse.cpp backendpetsc.cpp
  graphcsr.cpp  solvereigen.cpp
  matrixeigendense.cpp backendeigen.cpp vectoreigen.cpp matrixeigensparse.cpp
  solvernonlinear.cpp matrixblock.cpp vectorblock.cpp
  ${FEELPP_PETSC_SRCS}
  ${FEELPP_SLEPC_SRCS}
  ${FEELPP_TRILINOS_SRCS}
  )


file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
feel_append_src(feelalg "${SRC};${HDR}")

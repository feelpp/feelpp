/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2011-08-24

  Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file functionspetsc.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-08-24
 */
#ifndef __FunctionsPetsc_H
#define __FunctionsPetsc_H 1

#include <feel/feelcore/feel.hpp>

// Petsc include files.
#if defined( FEELPP_HAS_PETSC_H )

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
}
#else
# include <petscversion.h>
# include <petsc.h>
# include <petscsnes.h>
#endif

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/datamap.hpp>

namespace Feel
{

    //MatSolverPackageType matSolverPackageEnumType(std::string const& type );
void PetscPCFactorSetMatSolverPackage( PC & pc, MatSolverPackageType mspackt );

std::string PetscConvertKSPReasonToString( KSPConvergedReason reason );
std::string PetscConvertSNESReasonToString( SNESConvergedReason reason );

MatStructure PetscGetMatStructureEnum( Feel::MatrixStructure matStruc );

void PetscConvertIndexSplit( std::vector<IS> & isPetsc ,IndexSplit const& is,WorldComm const& worldcomm );
} // namespace Feel

#endif

#endif

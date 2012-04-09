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
   \file functionspetsc.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-08-24
 */

#include <feel/feelalg/functionspetsc.hpp>

namespace Feel
{

void
PetscPCFactorSetMatSolverPackage( PC & pc, MatSolverPackageType mspackt )
{
    int ierr = 0;
    Debug( 7010 ) << "[PetscPCFactorSetMatSolverPackage] :  " << mspackt << "\n";

    switch ( mspackt )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)

    case MATSOLVER_SPOOLES :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSPOOLES );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_SUPERLU :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSUPERLU );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_SUPERLU_DIST :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERSUPERLU_DIST );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_UMFPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERUMFPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_ESSL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERESSL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_LUSOL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERLUSOL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_MUMPS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERMUMPS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PASTIX :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPASTIX );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_MATLAB :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERMATLAB );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PETSC :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPETSC );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PLAPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERPLAPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 1)

    case MATSOLVER_BAS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MATSOLVERBAS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;
#endif

#else // PETSC < 3.2

    case MATSOLVER_SPOOLES :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SPOOLES );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_SUPERLU :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SUPERLU );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_SUPERLU_DIST :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_SUPERLU_DIST );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_UMFPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_UMFPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_ESSL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_ESSL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_LUSOL :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_LUSOL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_MUMPS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_MUMPS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PASTIX :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PASTIX );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_DSCPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_DSCPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_MATLAB :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_MATLAB );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PETSC :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PETSC );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case MATSOLVER_PLAPACK :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_PLAPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 1)

    case MATSOLVER_BAS :
        ierr = PCFactorSetMatSolverPackage( pc, ( char* ) MAT_SOLVER_BAS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;
#endif
#endif

    default:
        std::cerr << "ERROR:  Unsupported PETSC mat solver package: "
                  << mspackt               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

} // setPetscMatSolverPackageType


} // namespace Feel

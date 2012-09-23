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
#include <feel/feelcore/feelpetsc.hpp>
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

} // PetscPCFactorSetMatSolverPackage

std::string
PetscConvertKSPReasonToString( KSPConvergedReason reason )
{
    switch ( reason )
    {
        /* converged */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 2, 0 )
    case KSP_CONVERGED_RTOL_NORMAL     : return "CONVERGED_RTOL_NORMAL";
#endif
    case KSP_CONVERGED_RTOL            : return "CONVERGED_RTOL";
    case KSP_CONVERGED_ATOL            : return "CONVERGED_ATOL";
    case KSP_CONVERGED_ITS             : return "CONVERGED_ITS";
    case KSP_CONVERGED_CG_NEG_CURVE    : return "CONVERGED_CG_NEG_CURVE";
    case KSP_CONVERGED_CG_CONSTRAINED  : return "CONVERGED_CG_CONSTRAINED";
    case KSP_CONVERGED_STEP_LENGTH     : return "CONVERGED_STEP_LENGTH";
    case KSP_CONVERGED_HAPPY_BREAKDOWN : return "CONVERGED_HAPPY_BREAKDOWN";
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 2, 0 )
    case KSP_CONVERGED_ATOL_NORMAL     : return "CONVERGED_ATOL_NORMAL";
#endif

        /* diverged */
    case KSP_DIVERGED_NULL           : return "DIVERGED_NULL";
    case KSP_DIVERGED_ITS            : return "DIVERGED_ITS";
    case KSP_DIVERGED_DTOL           : return "DIVERGED_DTOL";
    case KSP_DIVERGED_BREAKDOWN      : return "DIVERGED_BREAKDOWN";
    case KSP_DIVERGED_BREAKDOWN_BICG : return "DIVERGED_BREAKDOWN_BICG";
    case KSP_DIVERGED_NONSYMMETRIC   : return "DIVERGED_NONSYMMETRIC";
    case KSP_DIVERGED_INDEFINITE_PC  : return "DIVERGED_INDEFINITE_PC";
    case KSP_DIVERGED_NAN            : return "DIVERGED_NAN";
    case KSP_DIVERGED_INDEFINITE_MAT : return "DIVERGED_INDEFINITE_MAT";

    case KSP_CONVERGED_ITERATING : return "CONVERGED_ITERATING";

    default: return "INDEFINE_KSP_REASON";

    }
} // PetscPCFactorSetMatSolverPackage

std::string
PetscConvertSNESReasonToString( SNESConvergedReason reason )
{
   switch ( reason )
    {
        /* converged */
    case SNES_CONVERGED_FNORM_ABS      : return "CONVERGED_FNORM_ABS";     // =  2, /* ||F|| < atol */
    case SNES_CONVERGED_FNORM_RELATIVE : return "CONVERGED_FNORM_RELATIVE";// =  3, /* ||F|| < rtol*||F_initial|| */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 3, 0 )
    case SNES_CONVERGED_SNORM_RELATIVE : return "CONVERGED_SNORM_RELATIVE";// =  4, /* Newton computed step size small; || delta x || < stol */
#endif
    case SNES_CONVERGED_ITS            : return "CONVERGED_ITS";           // =  5, /* maximum iterations reached */
    case SNES_CONVERGED_TR_DELTA       : return "CONVERGED_TR_DELTA";      // =  7,
        /* diverged */
    case SNES_DIVERGED_FUNCTION_DOMAIN : return "DIVERGED_FUNCTION_DOMAIN";// = -1, /* the new x location passed the function is not in the domain of F */
    case SNES_DIVERGED_FUNCTION_COUNT  : return "DIVERGED_FUNCTION_COUNT"; // = -2,
    case SNES_DIVERGED_LINEAR_SOLVE    : return "DIVERGED_LINEAR_SOLVE";   // = -3, /* the linear solve failed */
    case SNES_DIVERGED_FNORM_NAN       : return "DIVERGED_FNORM_NAN";      // = -4,
    case SNES_DIVERGED_MAX_IT          : return "DIVERGED_MAX_IT";         // = -5,
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 2, 0 )
    case SNES_DIVERGED_LINE_SEARCH     : return "DIVERGED_LINE_SEARCH";    // = -6, /* the line search failed */
#endif
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 3, 0 )
        case SNES_DIVERGED_INNER           : return "DIVERGED_INNER";          // = -7, /* inner solve failed */
#endif
    case SNES_DIVERGED_LOCAL_MIN       : return "DIVERGED_LOCAL_MIN";      // = -8, /* || J^T b || is small, implies converged to local minimum of F() */

    case SNES_CONVERGED_ITERATING      : return "CONVERGED_ITERATING";     // =  0

    default: return "INDEFINE_SNES_REASON";

    }
} //PetscConvertSNESReasonToString


MatStructure
PetscGetMatStructureEnum( Feel::MatrixStructure matStruc )
{
    switch (matStruc)
    {
    case Feel::SAME_NONZERO_PATTERN : return MatStructure::SAME_NONZERO_PATTERN;
    case Feel::DIFFERENT_NONZERO_PATTERN : return MatStructure::DIFFERENT_NONZERO_PATTERN;
    case Feel::SAME_PRECONDITIONER : return MatStructure::SAME_PRECONDITIONER;
    case Feel::SUBSET_NONZERO_PATTERN : return MatStructure::SUBSET_NONZERO_PATTERN;
        //case Feel::INVALID_STRUCTURE :
    default : return MatStructure::DIFFERENT_NONZERO_PATTERN;
    }
}

} // namespace Feel

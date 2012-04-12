/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-01-16

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file preconditionerpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-01-16
 */
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>

namespace Feel
{
template <typename T>
void PreconditionerPetsc<T>::apply( const Vector<T> & x, Vector<T> & y )
{
    VectorPetsc<T> & x_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( x ) );
    VectorPetsc<T> & y_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( y ) );

    Vec x_vec = x_pvec.vec();
    Vec y_vec = y_pvec.vec();

    int ierr = PCApply( M_pc,x_vec,y_vec );
    CHKERRABORT( this->M_comm,ierr );
}




template <typename T>
void PreconditionerPetsc<T>::init ()
{
    if ( !this->M_matrix )
    {
        std::cerr << "ERROR: No matrix set for PreconditionerPetsc, but init() called" << std::endl;
    }

    this->M_matrix->close();

    // Clear the preconditioner in case it has been created in the past
    if ( !this->M_is_initialized )
    {
        // Create the preconditioning object
        int ierr = PCCreate( this->M_comm,&M_pc );
        CHKERRABORT( this->M_comm,ierr );

        MatrixPetsc<T> * pmatrix = dynamic_cast<MatrixPetsc<T>*>( this->M_matrix.get() );

        M_mat = pmatrix->mat();
    }

    int ierr = PCSetOperators( M_pc,M_mat,M_mat,( MatStructure )SAME_NONZERO_PATTERN );
    CHKERRABORT( this->M_comm,ierr );

    // Set the PCType.  Note: this used to be done *before* the call to
    // PCSetOperators(), and only when !M_is_initialized, but
    // 1.) Some preconditioners (those employing sub-preconditioners,
    // for example) have to call PCSetUp(), and can only do this after
    // the operators have been set.
    // 2.) It should be safe to call set_petsc_preconditioner_type()
    // multiple times.
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc );

    this->M_is_initialized = true;
}


template <typename T>
void PreconditionerPetsc<T>::clear ()
{
    if ( this-> M_is_initialized )
    {
        this->M_is_initialized = false;
        PETSc::PCDestroy( M_pc );
    }
}




template <typename T>
void PreconditionerPetsc<T>::setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                                          const MatSolverPackageType & matSolverPackage_type,
                                                          PC & pc )
{
    mpi::communicator world;
    int ierr = 0;

    switch ( preconditioner_type )
    {
    case IDENTITY_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCNONE );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case CHOLESKY_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCCHOLESKY );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case ICC_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCICC );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case ILU_PRECOND:
    {
        // In serial, just set the ILU preconditioner type
        //if (Feel::n_processors() == 1)
        // change in parallel version
        if ( world.size() == 1 )
        {
            ierr = PCSetType ( pc, ( char* ) PCILU );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }

        else
        {
            // But PETSc has no truly parallel ILU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign ILU sub-preconditioners.
            ierr = PCSetType ( pc, ( char* ) PCBJACOBI );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );

            // Set ILU as the sub preconditioner type
            setPetscSubpreconditionerType( PCILU, pc );
        }

        break;
    }

    case LU_PRECOND:
    {

        // In serial, just set the LU preconditioner type
        //if (Feel::n_processors() == 1)
        // do be changed in parallel
        if ( world.size() == 1 || matSolverPackage_type == MATSOLVER_MUMPS)
        {
            ierr = PCSetType ( pc, ( char* ) PCLU );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }

        else
        {
            // But PETSc has no truly parallel LU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign LU sub-preconditioners.
            ierr = PCSetType ( pc, ( char* ) PCBJACOBI );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );

            // Set ILU as the sub preconditioner type
            setPetscSubpreconditionerType( PCLU, pc );
        }

        break;
    }

    case ASM_PRECOND:
    {
        // In parallel, I think ASM uses ILU by default as the sub-preconditioner...
        // I tried setting a different sub-preconditioner here, but apparently the matrix
        // is not in the correct state (at this point) to call PCSetUp().
        ierr = PCSetType ( pc, ( char* ) PCASM );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
    }

    case JACOBI_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCJACOBI );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCBJACOBI );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case SOR_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCSOR );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case EISENSTAT_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCEISENSTAT );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case AMG_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCHYPRE );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

#if !(PETSC_VERSION_LESS_THAN(2,1,2))

        // Only available for PETSC >= 2.1.2
    case USER_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCMAT );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;
#endif

    case SHELL_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCSHELL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case FIELDSPLIT_PRECOND: {
        ierr = PCSetType( pc,( char* ) PCFIELDSPLIT );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        const PCType subpctypes[5] = { PCLU,PCNONE,PCLU,PCLU,PCLU };
        const KSPType subksptypes[5] = { KSPPREONLY,KSPMINRES,KSPMINRES,KSPMINRES,KSPMINRES };
        setPetscFieldSplitPreconditionerType( PC_COMPOSITE_SCHUR, subksptypes, subpctypes, pc );
        break; }

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << preconditioner_type       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

    // Set additional options if we are doing AMG and
    // HYPRE is available
#ifdef FEELPP_HAS_PETSC_HYPRE

    if ( preconditioner_type == AMG_PRECOND )
    {
        ierr = PCHYPRESetType( pc, "boomeramg" );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

#endif

    // Let the commandline override stuff
    if ( preconditioner_type != AMG_PRECOND )
    {
        ierr = PCSetFromOptions( pc );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }
}



template <typename T>
#if PETSC_VERSION_LESS_THAN(3,0,0)
void PreconditionerPetsc<T>::setPetscSubpreconditionerType( PCType type, PC& pc )
#else
void PreconditionerPetsc<T>::setPetscSubpreconditionerType( const PCType type, PC& pc )
#endif
{
    // For catching PETSc error return codes
    int ierr = 0;

    // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
    // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
    //
    // "Object is in wrong state!"
    // "Matrix must be set first."
    //
    // error messages...
    ierr = PCSetUp( pc );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // To store array of local KSP contexts on this processor
    KSP* subksps;

    // the number of blocks on this processor
    int n_local;

    // The global number of the first block on this processor.
    // This is not used, so we just pass PETSC_NULL instead.
    // int first_local;

    // Fill array of local KSP contexts
    ierr = PCBJacobiGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Loop over sub-ksp objects, set ILU preconditioner
    for ( int i=0; i<n_local; ++i )
    {
        // Get pointer to sub KSP object's PC
        PC subpc;

        ierr = KSPGetPC( subksps[i], &subpc );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        // Set requested type on the sub PC
        ierr = PCSetType( subpc, type );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

}


template <typename T>
void
PreconditionerPetsc<T>::setPetscFieldSplitPreconditionerType( const PCCompositeType type,
                                                              const KSPType * subksptypes,
                                                              const PCType * subpctypes,
                                                              PC& pc )
{
    // For catching PETSc error return codes
    int ierr = 0;

    ierr = PCFieldSplitSetType( pc, type );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // call necessary before PCFieldSplitGetSubKSP
    ierr = PCSetUp( pc );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // the number of blocks on this processor
    int n_local;
    // To store array of local KSP contexts on this processor
    KSP* subksps;
    ierr = PCFieldSplitGetSubKSP(pc,&n_local,&subksps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Loop over sub-ksp objects, set ILU preconditioner
    for ( int i=0; i<n_local; ++i )
    {
        // Set requested type on the sub KSP
        ierr = KSPSetType ( subksps[i], subksptypes[i] );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        // Get pointer to sub KSP object's PC
        PC subpc;
        ierr = KSPGetPC( subksps[i], &subpc );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        // Set requested type on the sub PC
        ierr = PCSetType( subpc, subpctypes[i] );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

}




template class PreconditionerPetsc<double>;

}


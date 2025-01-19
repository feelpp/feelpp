
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008, 2009 Université de Grenoble 1

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
   \file solverlinearpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-27
 */
// $Id: petsc_linear_solver.C,v 1.5 2005/05/11 23:12:00 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include <feel/feelcore/feel.hpp>
#include <fmt/chrono.h>
#ifdef FEELPP_HAS_PETSC_H

#include <feel/feelalg/solverlinearpetsc.hpp>
#include <feel/feelalg/functionspetsc.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>


namespace Feel
{
extern "C"
{
    
#if PETSC_VERSION_LESS_THAN(2,2,1)
    typedef int PetscErrorCode;
    typedef int PetscInt;
#endif
    PetscErrorCode __feel_petsc_monitor(KSP ksp,PetscInt it,PetscReal rnorm,void* ctx)
    {
        SolverLinear<double> *s  = static_cast<SolverLinear<double>*>( ctx );
        if ( !s ) return 1;
        if ( s->worldComm().isMasterRank() )
            std::cout << fmt::format( "[{:%Y-%m-%d %H:%M:%S} - [{}] ] #{} KSP Residual norm {:.4e}", fmt::localtime( std::time(nullptr) ), s->prefix(), it, rnorm ) << std::endl;
//        if ( s->worldComm().isMasterRank() )
//            std::cout << " " << it  << " " << s->prefix() << " KSP Residual norm " << std::scientific << rnorm << "\n";
        return 0;
    }
#if PETSC_VERSION_LESS_THAN(3,0,1)
    PetscErrorCode __feel_petsc_preconditioner_setup ( void * ctx )
    {
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );
        preconditioner->init();

        VLOG(2) << "__feel_petsc_preconditioner_setup:: init prec\n";

        return 0;
    }

    PetscErrorCode __feel_petsc_preconditioner_apply( void *ctx, Vec x, Vec y )
    {
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );

        VectorPetsc<double> x_vec( x );
        VectorPetsc<double> y_vec( y );

        preconditioner->apply( x_vec,y_vec );

        return 0;
    }
    PetscErrorCode __feel_petsc_preconditioner_view( void *ctx, PetscViewer viewer)
    {
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );
        preconditioner->view();
        return 0;
    }

#else
    PetscErrorCode __feel_petsc_preconditioner_setup ( PC pc )
    {
        void *ctx;
        PetscErrorCode ierr = PCShellGetContext( pc,&ctx );
        CHKERRQ( ierr );
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
        bool reusePrec = preconditioner->reusePrec();
        // if we are here and reusePrec option, need to rebuild the preconditioner
        if ( reusePrec )
            preconditioner->setPrecMatrixStructure( MatrixStructure::SAME_NONZERO_PATTERN );
#endif
        // build preconditioner
        preconditioner->init();

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
        // tell to not rebuild the preconditioner after
        if ( reusePrec )
            preconditioner->setPrecMatrixStructure( MatrixStructure::SAME_PRECONDITIONER );
#endif
        VLOG(2) << "__feel_petsc_preconditioner_setup: init prec " << preconditioner->name() << "\n";
        return 0;
    }

    PetscErrorCode __feel_petsc_preconditioner_apply( PC pc, Vec x, Vec y )
    {
        void *ctx;
        PetscErrorCode ierr = PCShellGetContext( pc,&ctx );
        CHKERRQ( ierr );
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );

        PreconditionerPetsc<double> * preconditionerPetsc = dynamic_cast<PreconditionerPetsc<double>*>( preconditioner );
        if ( preconditionerPetsc != NULL )
        {
            preconditionerPetsc->apply(x,y);
        }
        else
        {
            std::shared_ptr<VectorPetsc<double> > x_vec;
            std::shared_ptr<VectorPetsc<double> > y_vec;
            if ( preconditioner->worldComm().localSize() > 1 )
            {
                CHECK ( preconditioner->matrix() ) << "matrix is not defined";
                Vec lx, ly;
                VecGhostGetLocalForm(x,&lx);
                VecGhostGetLocalForm(y,&ly);
                if ( lx )
                    x_vec.reset( new VectorPetscMPI<double>( x, preconditioner->matrix()->mapColPtr() ) );
                else
                    x_vec.reset( new VectorPetscMPIRange<double>( x, preconditioner->matrix()->mapColPtr() ) );
                if ( ly )
                    y_vec.reset( new VectorPetscMPI<double>( y, preconditioner->matrix()->mapRowPtr() ) );
                else
                    y_vec.reset( new VectorPetscMPIRange<double>( y, preconditioner->matrix()->mapRowPtr() ) );
                VecGhostRestoreLocalForm(x,&lx);
                VecGhostRestoreLocalForm(y,&ly);
            }
            else
            {
                x_vec.reset( new VectorPetsc<double>( x,preconditioner->matrix()->mapColPtr() ) );
                y_vec.reset( new VectorPetsc<double>( y,preconditioner->matrix()->mapRowPtr() ) );
            }
            preconditioner->apply( *x_vec,*y_vec );
        }

        return 0;
    }
    PetscErrorCode __feel_petsc_preconditioner_view( PC pc, PetscViewer viewer)
    {
        void *ctx;
        PetscErrorCode ierr = PCShellGetContext( pc,&ctx );
        CHKERRQ( ierr );
        Preconditioner<double> * preconditioner = static_cast<Preconditioner<double>*>( ctx );
        preconditioner->view();
        return 0;
    }
#endif
} // end extern "C"

/*----------------------- functions ----------------------------------*/
template <typename T>
void
SolverLinearPetsc<T>::clear ()
{
    PetscBool pinit;
    PetscInitialized( &pinit );
    if ( pinit && this->initialized() )
    {
        this->setInitialized( false );

        int ierr=0;

        // 2.1.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

        ierr = SLESDestroy( M_sles );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // 2.2.0 & newer style
#else
        FEELPP_ASSERT( M_ksp != 0 ).error( "invalid ksp" );
        ierr = PETSc::KSPDestroy( M_ksp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
#endif

        // Mimic PETSc default solver and preconditioner
        this->setSolverType(  GMRES );

        if ( this->worldComm().globalComm().size() == 1 )
            this->setPreconditionerType( LU_PRECOND );

        else
            this->setPreconditionerType( BLOCK_JACOBI_PRECOND );
    }
}



template <typename T>
void SolverLinearPetsc<T>::init ()
{
    // Initialize the data structures if not done so already.
    if ( !this->initialized() )
    {
        this->setInitialized(  true );

        int ierr=0;

        // 2.1.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

        // Create the linear solver context
        ierr = SLESCreate ( this->worldComm().globalComm(), &M_sles );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Create the Krylov subspace & preconditioner contexts
        ierr = SLESGetKSP       ( M_sles, &M_ksp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        ierr = SLESGetPC        ( M_sles, &M_pc );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Have the Krylov subspace method use our good initial guess rather than 0
        ierr = KSPSetInitialGuessNonzero ( M_ksp, PETSC_TRUE );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Set user-specified  solver and preconditioner types
        this->setPetscSolverType();
        this->setPetscPreconditionerType();
        this->setPetscConstantNullSpace();

        // Set the options from user-input
        // Set runtime options, e.g.,
        //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        //  These options will override those specified above as long as
        //  SLESSetFromOptions() is called _after_ any other customization
        //  routines.

        ierr = SLESSetFromOptions ( M_sles );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // 2.2.0 & newer style
#else

        // Create the linear solver context
        ierr = KSPCreate ( this->worldComm().globalComm(), &M_ksp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Create the preconditioner context
        ierr = KSPGetPC        ( M_ksp, &M_pc );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Have the Krylov subspace method use our good initial guess rather than 0
        ierr = KSPSetInitialGuessNonzero ( M_ksp, this->M_kspUseInitialGuessNonZero?PETSC_TRUE:PETSC_FALSE );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        // Set user-specified  solver and preconditioner types
        this->setPetscSolverType();
        this->setPetscConstantNullSpace();


        // Set the options from user-input
        // Set runtime options, e.g.,
        //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        //  These options will override those specified above as long as
        //  KSPSetFromOptions() is called _after_ any other customization
        //  routines.
        //ierr = PCSetFromOptions ( M_pc );
        //CHKERRABORT( this->worldComm().globalComm(),ierr );
        ierr = KSPSetFromOptions ( M_ksp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );


#endif

        // Have the Krylov subspace method use our good initial guess
        // rather than 0, unless the user requested a KSPType of
        // preonly, which complains if asked to use initial guesses.
#if PETSC_VERSION_LESS_THAN(3,0,0)
        KSPType ksp_type;
#else
#if PETSC_VERSION_LESS_THAN(3,4,0)
        const KSPType ksp_type;
#else
        KSPType ksp_type;
#endif
#endif

        ierr = KSPGetType ( M_ksp, &ksp_type );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        if ( std::string((char*)ksp_type) == std::string( ( char* )KSPPREONLY ) )
        {
            ierr = KSPSetInitialGuessNonzero ( M_ksp, PETSC_FALSE );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
        }
        else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPGMRES ) )
        {
            ierr = KSPGMRESSetRestart( M_ksp, this->M_nRestartGMRES );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
        }
        else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPFGMRES ) )
        {
            ierr = KSPGMRESSetRestart( M_ksp, this->M_nRestartFGMRES );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
            if ( this->M_preconditioner )
                this->M_preconditioner->setSide( preconditioner_type::RIGHT );
        }
        else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPGCR ) )
        {
            ierr = KSPGCRSetRestart( M_ksp, this->M_nRestartGCR );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
            if ( this->M_preconditioner )
                this->M_preconditioner->setSide( preconditioner_type::RIGHT );
        }
        // Notify PETSc of location to store residual history.
        // This needs to be called before any solves, since
        // it sets the residual history length to zero.  The default
        // behavior is for PETSc to allocate (internally) an array
        // of size 1000 to hold the residual norm history.
        ierr = KSPSetResidualHistory( M_ksp,
                                      PETSC_IGNORE,   // pointer to the array which holds the history
                                      PETSC_DECIDE, // size of the array holding the history
                                      PETSC_TRUE ); // Whether or not to reset the history for each solve.
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        //If there is a preconditioner object we need to set the internal setup and apply routines
        if ( this->M_preconditioner )
        {
            VLOG(2) << "preconditioner: "  << this->M_preconditioner << "\n";

            PCSetType(M_pc, PCSHELL);
            PCShellSetName( M_pc, this->M_preconditioner->name().c_str() );
            PCShellSetContext( M_pc,( void* )this->M_preconditioner.get() );
            PCShellSetSetUp( M_pc,__feel_petsc_preconditioner_setup );
            PCShellSetApply( M_pc,__feel_petsc_preconditioner_apply );
            PCShellSetView( M_pc,__feel_petsc_preconditioner_view );
#if PETSC_VERSION_LESS_THAN(3,2,0)
            const PCType pc_type;
#else
            PCType pc_type;
#endif
            ierr = PCGetType ( M_pc, &pc_type );
            CHKERRABORT( this->worldComm().globalComm(),ierr );

            switch( this->M_preconditioner->side() )
            {
            default:
            case preconditioner_type::LEFT:
                VLOG(2) << " . PC is set to left side\n";
#if PETSC_VERSION_LESS_THAN(3,2,0)
                KSPSetPreconditionerSide( M_ksp, PC_LEFT );
#else
                KSPSetPCSide( M_ksp, PC_LEFT );
#endif
                break;
            case preconditioner_type::RIGHT:
                VLOG(2) << " . PC is set to right side\n";
#if PETSC_VERSION_LESS_THAN(3,2,0)
                KSPSetPreconditionerSide( M_ksp, PC_RIGHT );
#else
                KSPSetPCSide( M_ksp, PC_RIGHT );
#endif
                break;
            case preconditioner_type::SYMMETRIC:
                VLOG(2) << " . PC is set to symmetric\n";
#if PETSC_VERSION_LESS_THAN(3,4,0)
                KSPSetPreconditionerSide( M_ksp, PC_SYMMETRIC );
#else
                KSPSetPCSide( M_ksp, PC_SYMMETRIC );
#endif
                break;
            }

            VLOG(2) << "preconditioner set as "  << pc_type << "\n";
        }
        else
        {
            this->setPetscPreconditionerType();
            // sets the software that is used to perform the factorization
            PetscPCFactorSetMatSolverPackage( M_pc,this->matSolverPackageType() );
        }

        if ( this->M_showKSPMonitor )
        {
            //KSPMonitorSet( M_ksp,KSPMonitorDefault,PETSC_IGNORE,PETSC_IGNORE );
            KSPMonitorSet( M_ksp,__feel_petsc_monitor,(void*) this,PETSC_IGNORE );
        }

        // The value can be checked with --(prefix.)ksp-view=1
        this->check( KSPSetNormType(M_ksp, kspNormTypeConvertStrToEnum( this->M_kspNormType ) ) );

    }
}








template <typename T>
typename SolverLinearPetsc<T>::solve_return_type
SolverLinearPetsc<T>::solve ( MatrixSparse<T> const&  matrix_in,
                              MatrixSparse<T> const&  precond_in,
                              Vector<T> & solution_in,
                              Vector<T> const& rhs_in,
                              const double tol,
                              const unsigned int m_its,
                              bool transpose )
{
    this->setWorldComm( matrix_in.worldCommPtr() );
    this->init ();

    MatrixPetsc<T> * matrix   = const_cast<MatrixPetsc<T> *>( dynamic_cast<MatrixPetsc<T> const*>( &matrix_in ) );
    MatrixPetsc<T> * precond  = const_cast<MatrixPetsc<T> *>( dynamic_cast<MatrixPetsc<T> const*>( &precond_in ) );
    VectorPetsc<T> * solution = dynamic_cast<VectorPetsc<T>*>( &solution_in );
    VectorPetsc<T> * rhs      = const_cast<VectorPetsc<T> *>( dynamic_cast<VectorPetsc<T> const*>( &rhs_in ) );

    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
    DCHECK( matrix   != nullptr ) << "non petsc matrix structure";
    DCHECK( precond  != nullptr ) << "non petsc matrix structure";
    DCHECK( solution != nullptr ) << "non petsc vector structure";
    DCHECK( rhs      != nullptr ) << "non petsc vector structure";

    int ierr=0;
    int its=0;
    PetscReal final_resid=0.;

    solution->close ();


    if ( !this->M_preconditioner && this->preconditionerType() == FIELDSPLIT_PRECOND )
        matrix->updatePCFieldSplit( M_pc );

    if ( this->M_nullSpace && this->M_nullSpace->size() > 0 )
        this->updateNullSpace( matrix->mat(), rhs->vec() );
    if ( this->M_nearNullSpace && this->M_nearNullSpace->size() > 0 )
        this->updateNearNullSpace( matrix->mat() );

    //   // If matrix != precond, then this means we have specified a
    //   // special preconditioner, so reset preconditioner type to PCMAT.
    //   if (matrix != precond)
    //     {
    //       this->_preconditioner_type = USER_PRECOND;
    //       this->set_petsc_preconditioner_type ();
    //     }

    
    // 2.1.x & earlier style
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

    // Set operators. The input matrix works as the preconditioning matrix
    ierr = SLESSetOperators( M_sles, matrix->mat(), precond->mat(),
                             SAME_NONZERO_PATTERN );
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values.
    ierr = KSPSetTolerances ( M_ksp,
                              this->rTolerance(),
                              this->aTolerance(),
                              this->dTolerance(),
                              this->maxIterations() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // makes the default convergence test use || B*(b - A*(initial guess))||
    // instead of || B*b ||. In the case of right preconditioner or if
    // KSPSetNormType(ksp,KSP_NORM_UNPRECONDIITONED) is used there is no B in
    // the above formula. UIRNorm is short for Use Initial Residual Norm.
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,4,4)
    KSPConvergedDefaultSetUIRNorm( M_ksp );
#else
    KSPDefaultConvergedSetUIRNorm( M_ksp );
#endif


    // Solve the linear system
    ierr = SLESSolve ( M_sles, rhs->vec(), solution->vec(), &its );
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm ( M_ksp, &final_resid );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.2.0
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 2) && (PETSC_VERSION_SUBMINOR == 0)

    // Set operators. The input matrix works as the preconditioning matrix
    ierr = KSPSetOperators( M_ksp, matrix->mat(), precond->mat(),
                            MatStructure::SAME_NONZERO_PATTERN );
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values.
    // Convergence is detected at iteration k if
    // ||r_k||_2 < max(rtol*||b||_2 , abstol)
    // where r_k is the residual vector and b is the right-hand side.  Note that
    // it is the *maximum* of the two values, the larger of which will almost
    // always be rtol*||b||_2.
    ierr = KSPSetTolerances ( M_ksp,
                              this->rTolerance(),
                              this->aTolerance(),
                              this->dTolerance(),
                              this->maxIterations() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    // Set the solution vector to use
    ierr = KSPSetSolution ( M_ksp, solution->vec() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Set the RHS vector to use
    ierr = KSPSetRhs ( M_ksp, rhs->vec() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // makes the default convergence test use || B*(b - A*(initial guess))||
    // instead of || B*b ||. In the case of right preconditioner or if
    // KSPSetNormType(ksp,KSP_NORM_UNPRECONDIITONED) is used there is no B in
    // the above formula. UIRNorm is short for Use Initial Residual Norm.
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,4,4)
    KSPConvergedDefaultSetUIRNorm( M_ksp );
#else
    KSPDefaultConvergedSetUIRNorm( M_ksp );
#endif

    // Solve the linear system
    if ( transpose )
        ierr = KSPSolveTranspose ( M_ksp );

    else
        ierr = KSPSolve ( M_ksp );

    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber ( M_ksp, &its );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm ( M_ksp, &final_resid );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.2.1 & newer style
#else
    //std::cout << "sles: " << this->precMatrixStructure() << "\n";
    // Set operators. The input matrix works as the preconditioning matrix
#if PETSC_VERSION_LESS_THAN(3,5,0)
    ierr = KSPSetOperators( M_ksp, matrix->mat(), precond->mat(),
                            PetscGetMatStructureEnum(this->precMatrixStructure()) );
#else
    ierr = KSPSetReusePreconditioner( M_ksp, (this->precMatrixStructure() == Feel::SAME_PRECONDITIONER)? PETSC_TRUE : PETSC_FALSE );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    ierr = KSPSetOperators( M_ksp, matrix->mat(), precond->mat() );
#endif
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values.
    ierr = KSPSetTolerances ( M_ksp,
                              this->rTolerance(),
                              //1e-15,
                              this->aTolerance(),
                              this->dTolerance(),
                              this->maxIterations() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    //PreconditionerPetsc<T>::setPetscPreconditionerType( this->preconditionerType(),this->matSolverPackageType(),M_pc, this->worldComm() );


    // makes the default convergence test use || B*(b - A*(initial guess))||
    // instead of || B*b ||. In the case of right preconditioner or if
    // KSPSetNormType(ksp,KSP_NORM_UNPRECONDIITONED) is used there is no B in
    // the above formula. UIRNorm is short for Use Initial Residual Norm.
#if PETSC_VERSION_LESS_THAN(3,5,0)
    KSPDefaultConvergedSetUIRNorm( M_ksp );
#else
    KSPConvergedDefaultSetUIRNorm( M_ksp );
#endif

    // Solve the linear system
    if ( transpose )
        ierr = KSPSolveTranspose ( M_ksp, rhs->vec(), solution->vec() );

    else
        ierr = KSPSolve ( M_ksp, rhs->vec(), solution->vec() );

    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Get the number of iterations required for convergence
    ierr = KSPGetIterationNumber ( M_ksp, &its );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Get the norm of the final residual to return to the user.
    ierr = KSPGetResidualNorm ( M_ksp, &final_resid );
    //std::cout << "final residual = " << final_resid << "\n";
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    KSPConvergedReason reason;
    KSPGetConvergedReason( M_ksp,&reason );

    if ( this->M_kspView ) //boption( _prefix=this->prefix(), _name="ksp-view" ) )
        check( KSPView( M_ksp, PETSC_VIEWER_STDOUT_WORLD ) );

    if ( reason==KSP_DIVERGED_INDEFINITE_PC )
    {
        LOG(INFO) << "[solverlinearpetsc] Divergence because of indefinite preconditioner;\n";
        LOG(INFO) << "[solverlinearpetsc] Run the executable again but with '-pc_factor_shift_type POSITIVE_DEFINITE' option.\n";
    }

    else if ( reason<0 )
    {
        LOG(INFO) <<"[solverlinearpetsc] Other kind of divergence: this should not happen.\n";
    }

    bool hasConverged;

    if ( reason> 0 )
        {
            hasConverged=true;
            if (this->showKSPConvergedReason() && this->worldComm().globalRank() == this->worldComm().masterRank() )
                std::cout<< "Linear solve converged due to " << PetscConvertKSPReasonToString(reason)
                         << " iterations " << its << std::endl;
        }
    else
        {
            hasConverged=false;
            if (this->showKSPConvergedReason() && this->worldComm().globalRank() == this->worldComm().masterRank() )
                std::cout<< "Linear solve did not converge due to " << PetscConvertKSPReasonToString(reason)
                         << " iterations " << its << std::endl;
        }

#endif
    // return the # of its. and the final residual norm.
    //return std::make_pair(its, final_resid);
    return solve_return_type( boost::make_tuple( hasConverged, its, final_resid ) );


}

template <typename T>
typename SolverLinearPetsc<T>::solve_return_type
SolverLinearPetsc<T>::solve ( MatrixShell<T>  const &mat,
                              Vector<T> & x,
                              Vector<T> const& b,
                              const double tolerance,
                              const unsigned int maxit,
                              bool transpose )
{
    LOG(ERROR) << "invalid call to solve() using matshell\n";
    return solve_return_type( boost::make_tuple( false, 0, 0 ) );
}

template <typename T>
void
SolverLinearPetsc<T>::getResidualHistory( std::vector<double>& hist )
{
    int ierr = 0;
    int its  = 0;

    // Fill the residual history vector with the residual norms
    // Note that GetResidualHistory() does not copy any values, it
    // simply sets the pointer p.  Note that for some Krylov subspace
    // methods, the number of residuals returned in the history
    // vector may be different from what you are expecting.  For
    // example, TFQMR returns two residual values per iteration step.
#if PETSC_VERSION_LESS_THAN(3,15,0)
    double* p;
#else
    const double* p;
#endif
    ierr = KSPGetResidualHistory( M_ksp, &p, &its );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Check for early return
    if ( its == 0 ) return;

    // Create space to store the result
    hist.resize( its );

    // Copy history into the vector provided by the user.
    for ( int i=0; i<its; ++i )
    {
        hist[i] = *p;
        p++;
    }
}




template <typename T>
typename SolverLinearPetsc<T>::real_type
SolverLinearPetsc<T>::getInitialResidual()
{
    int ierr = 0;
    int its  = 0;

    // Fill the residual history vector with the residual norms
    // Note that GetResidualHistory() does not copy any values, it
    // simply sets the pointer p.  Note that for some Krylov subspace
    // methods, the number of residuals returned in the history
    // vector may be different from what you are expecting.  For
    // example, TFQMR returns two residual values per iteration step.
#if PETSC_VERSION_LESS_THAN(3,15,0)
    double* p;
#else
    const double* p;
#endif
    ierr = KSPGetResidualHistory( M_ksp, &p, &its );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Check no residual history
    if ( its == 0 )
    {
        std::cerr << "No iterations have been performed, returning 0." << std::endl;
        return 0.;
    }

    // Otherwise, return the value pointed to by p.
    return *p;
}


template <typename T>
void
SolverLinearPetsc<T>::setPetscConstantNullSpace()
{
    int ierr = 0;
    if ( this->M_nullSpace && this->M_nullSpace->size() > 0 )
    {
        //std::cout << "define nullspace in petsc with size "<< this->M_nullSpace->size() <<"\n";
    }
    else if ( M_constant_null_space )
    {
        std::cout << "use nullspace in petsc\n";
        MatNullSpace nullsp;

        ierr = MatNullSpaceCreate( PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_IGNORE, &nullsp );
        CHKERRABORT( this->worldComm().globalComm(), ierr );
#if PETSC_VERSION_LESS_THAN( 3,5,4 )
        ierr = KSPSetNullSpace( M_ksp, nullsp );
#else
        Mat A;
        ierr = KSPGetOperators( M_ksp, &A, NULL );
        CHKERRABORT( this->worldComm().globalComm(), ierr );
        ierr = MatSetNullSpace( A, nullsp );
#endif
        CHKERRABORT( this->worldComm().globalComm(), ierr );
        PETSc::MatNullSpaceDestroy( nullsp );
    }
}

template <typename T>
void
SolverLinearPetsc<T>::updateNullSpace( Mat A, Vec rhs )
{
    if ( !this->M_nullSpace ) return;
    if ( this->M_nullSpace->size() == 0 ) return;

    int ierr = 0;
    int dimNullSpace = this->M_nullSpace->size();
    std::vector<Vec> petsc_vec(dimNullSpace);
    for ( int k = 0 ; k<dimNullSpace ; ++k )
        petsc_vec[k] =  dynamic_cast<const VectorPetsc<T>*>( &this->M_nullSpace->basisVector(k) )->vec();

#if 0
    // reorthornomalisation with petsc
    PetscScalar dots[5];
    for (int i=0/*dim*/; i<dimNullSpace; i++) {
        /* Orthonormalize vec[i] against vec[0:i-1] */
        VecMDot(petsc_vec[i],i,petsc_vec.data(),dots);
        for (int j=0; j<i; j++) dots[j] *= -1.;
        VecMAXPY(petsc_vec[i],i,dots,petsc_vec.data()/*vec*/);
        VecNormalize(petsc_vec[i],NULL);
    }
#endif


#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    MatNullSpace nullsp;
    ierr = MatNullSpaceCreate( this->worldComm(), PETSC_FALSE , dimNullSpace, petsc_vec.data()/*PETSC_IGNORE*/, &nullsp );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    //ierr = MatNullSpaceView( nullsp, PETSC_VIEWER_STDOUT_WORLD );
    //CHKERRABORT( this->worldComm().globalComm(),ierr );

    ierr = MatSetNullSpace(A,nullsp);
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    //ierr = MatNullSpaceRemove(nullsp,rhs);
    //CHKERRABORT( this->worldComm().globalComm(),ierr );

    bool checkNullSpace = false;
    if ( checkNullSpace )
    {
        PetscBool isNull;
        ierr = MatNullSpaceTest(nullsp, A, &isNull);
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        CHECK( isNull ) << "nullspace is not apply on this matrix";
    }

    PETSc::MatNullSpaceDestroy( nullsp );
#endif
}

template <typename T>
void
SolverLinearPetsc<T>::updateNearNullSpace( Mat A )
{
    if ( !this->M_nearNullSpace ) return;
    if ( this->M_nearNullSpace->size() == 0 ) return;

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    int ierr = 0;
    int dimNullSpace = this->M_nearNullSpace->size();
    std::vector<Vec> petsc_vec(dimNullSpace);
    for ( int k = 0 ; k<dimNullSpace ; ++k )
        petsc_vec[k] =  dynamic_cast<const VectorPetsc<T>*>( &this->M_nearNullSpace->basisVector(k) )->vec();
    MatNullSpace nullsp;
    ierr = MatNullSpaceCreate( this->worldComm(), PETSC_FALSE , dimNullSpace, petsc_vec.data()/*PETSC_IGNORE*/, &nullsp );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    ierr = MatSetNearNullSpace( A, nullsp);
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    PETSc::MatNullSpaceDestroy( nullsp );
#endif
}
template <typename T>
void
SolverLinearPetsc<T>::setPetscSolverType()
{
    int ierr = 0;
    DVLOG(2) << "[SolverLinearPetsc] solver type: " << this->solverType() << "\n";

    switch (this->solverType())
    {
        case CG:
            ierr = KSPSetType(M_ksp, (char*)KSPCG);
            break;

        case CR:
            ierr = KSPSetType(M_ksp, (char*)KSPCR);
            break;

        case CGS:
            ierr = KSPSetType(M_ksp, (char*)KSPCGS);
            break;

        case BICG:
            ierr = KSPSetType(M_ksp, (char*)KSPBICG);
            break;

        case TCQMR:
            ierr = KSPSetType(M_ksp, (char*)KSPTCQMR);
            break;

        case TFQMR:
            ierr = KSPSetType(M_ksp, (char*)KSPTFQMR);
            break;

        case LSQR:
            ierr = KSPSetType(M_ksp, (char*)KSPLSQR);
            break;

        case BICGSTAB:
            ierr = KSPSetType(M_ksp, (char*)KSPBCGS);
            break;

        case MINRES:
            ierr = KSPSetType(M_ksp, (char*)KSPMINRES);
            break;

        case GMRES:
            ierr = KSPSetType(M_ksp, (char*)KSPGMRES);
            break;

        case FGMRES:
            ierr = KSPSetType(M_ksp, (char*)KSPFGMRES);
            break;

        case RICHARDSON:
            ierr = KSPSetType(M_ksp, (char*)KSPRICHARDSON);
            break;

        case CHEBYSHEV:
            ierr = KSPSetType(M_ksp, (char*)KSPCHEBYSHEV);
            break;

        case PREONLY:
            ierr = KSPSetType(M_ksp, (char*)KSPPREONLY);
            break;

        case GCR:
            ierr = KSPSetType(M_ksp, (char*)KSPGCR);
            break;

        // new KSP solvers
        case CGN:        // CG on Normal Equations => KSPCGNE
            ierr = KSPSetType(M_ksp, (char*)KSPCGNE);
            break;

        case QMR:        // NOTE: PETSc has "KSPQMRCGS" (qmrcgs) 
                         // There is no direct "qmr" in modern PETSc 
                         // Possibly fallback to KSPQMRCGS or 
                         // just ignore if not used
            ierr = KSPSetType(M_ksp, (char*)KSPQMRCGS);
            break;

        case JACOBI:      // Historically, Jacobi is not a KSP in PETSc, 
                          // but if we want to emulate, we might do:
            // We have to fake it or do a preonly+pc=jacobi
            ierr = KSPSetType(M_ksp, (char*)KSPPREONLY);
            break;

        case SOR_FORWARD: // Again, no direct KSP in PETSc
        case SOR_BACKWARD:
        case SSOR:
            // Typically we do KSPPREONLY + PCType SOR
            ierr = KSPSetType(M_ksp, (char*)KSPPREONLY);
            break;

        case DGMRES:
            ierr = KSPSetType(M_ksp, (char*)KSPDGMRES);
            break;

        case LGMRES:
            ierr = KSPSetType(M_ksp, (char*)KSPLGMRES);
            break;

        case PGMRES:
            ierr = KSPSetType(M_ksp, (char*)KSPPGMRES);
            break;

        case PIPEGMRES:
            // Could map to "KSPPGMRES" or "KSPPIPEFGMRES" 
            // depending on which pipeline GMRES is intended
            ierr = KSPSetType(M_ksp, (char*)KSPPIPEFGMRES);
            break;

        case PIPECG:
            ierr = KSPSetType(M_ksp, (char*)KSPPIPECG);
            break;

        case PIPECR:
            ierr = KSPSetType(M_ksp, (char*)KSPPIPECR);
            break;

        case BCGSL:
            ierr = KSPSetType(M_ksp, (char*)KSPBCGSL);
            break;

        case FBCGS:
            ierr = KSPSetType(M_ksp, (char*)KSPFBCGS);
            break;

        case IBCGS:
            ierr = KSPSetType(M_ksp, (char*)KSPIBCGS);
            break;

        case FCG:
            ierr = KSPSetType(M_ksp, (char*)KSPFCG);
            break;

        case FBCGSR:
            ierr = KSPSetType(M_ksp, (char*)KSPFBCGSR);
            break;

        case GROPPCG:
            ierr = KSPSetType(M_ksp, (char*)KSPGROPPCG);
            break;

        case PIPECGRR:
            ierr = KSPSetType(M_ksp, (char*)KSPPIPECGRR);
            break;

        case PIPEFCG:
            ierr = KSPSetType(M_ksp, (char*)KSPPIPEFCG);
            break;

        case CGLS:
            ierr = KSPSetType(M_ksp, (char*)KSPCGLS);
            break;

        case NASH:
            ierr = KSPSetType(M_ksp, (char*)KSPNASH);
            break;

        case STCG:
            ierr = KSPSetType(M_ksp, (char*)KSPSTCG);
            break;

        case GLTR:
            ierr = KSPSetType(M_ksp, (char*)KSPGLTR);
            break;

        case QCG:
            ierr = KSPSetType(M_ksp, (char*)KSPQCG);
            break;

        case FETIDP:
            ierr = KSPSetType(M_ksp, (char*)KSPFETIDP);
            break;

        case TSIRM:
            ierr = KSPSetType(M_ksp, (char*)KSPTSIRM);
            break;

        case SYMMLQ:
            ierr = KSPSetType(M_ksp, (char*)KSPSYMMLQ);
            break;

        case PYTHON:
            ierr = KSPSetType(M_ksp, (char*)KSPPYTHON);
            break;

        case NONE:
            ierr = KSPSetType(M_ksp, (char*)KSPNONE);
            break;

        case INVALID_SOLVER:
        default:
            std::cerr << "ERROR: Unsupported PETSC SolverType: " 
                      << this->solverType() << "\n"
                      << "Continuing with PETSc default (gmres)\n";
            ierr = KSPSetType(M_ksp, (char*)KSPGMRES);
            break;
    } // end switch

    CHKERRABORT(this->worldComm().globalComm(), ierr);
} // setPetscSolverType()






template <typename T>
void
SolverLinearPetsc<T>::setPetscPreconditionerType()
{
    int ierr = 0;
    DVLOG(2) << "[SolverLinearPetsc] preconditioner type: " << this->preconditionerType() << "\n";

    // Attempt default factor solver package if none given
    #if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,9,0)
    ierr = PCFactorSetMatSolverType( M_pc, MATSOLVERUMFPACK );
    if ( ierr ) ierr = PCFactorSetMatSolverType( M_pc, MATSOLVERSUPERLU );
    if ( ierr ) ierr = PCFactorSetMatSolverType( M_pc, MATSOLVERPETSC );
    #elif PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,2,0)
    ierr = PCFactorSetMatSolverPackage( M_pc, MATSOLVERUMFPACK );
    if ( ierr ) ierr = PCFactorSetMatSolverPackage( M_pc, MATSOLVERSUPERLU );
    if ( ierr ) ierr = PCFactorSetMatSolverPackage( M_pc, MATSOLVERPETSC );
    #endif

    switch ( this->preconditionerType() )
    {
      case IDENTITY_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCNONE );
          break;

      case CHOLESKY_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCCHOLESKY );
          break;

      case ICC_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCICC );
          break;

      case ILU_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCILU );
          break;

      case LU_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCLU );
          break;

      case ASM_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCASM );
          break;

      case GASM_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCGASM );
          break;

      case JACOBI_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCJACOBI );
          break;

      case BLOCK_JACOBI_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCBJACOBI );
          break;

      case SOR_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSOR );
          break;

      case SSOR_PRECOND:
          // PETSc typically does SOR with -pc_sor_symmetric 
          // or user sets if they want SSOR specifically
          ierr = PCSetType( M_pc, (char*)PCSOR );
          break;

      case EISENSTAT_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCEISENSTAT );
          break;

      case USER_PRECOND:
          // "USER_PRECOND" might have meant PCMAT historically
          ierr = PCSetType( M_pc, (char*)PCMAT );
          break;

      case SHELL_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSHELL );
          break;

      case FIELDSPLIT_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCFIELDSPLIT );
          // optionally use a default Schur factorization
          ierr = PCFieldSplitSetType( M_pc, PC_COMPOSITE_SCHUR );
          break;

      case ML_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCML );
          break;

      case GAMG_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCGAMG );
          break;

      case BOOMERAMG_PRECOND:
      case AMS_PRECOND:
          // Typically we do PCSetType( m_pc, PCHYPRE ), then set 
          // -pc_hypre_type boomeramg or -pc_hypre_type ams at runtime 
          ierr = PCSetType( M_pc, (char*)PCHYPRE );
          break;

      case REDUNDANT_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCREDUNDANT );
          break;

      case NONE_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCNONE );
          break;

      // New PC entries 

      case BDDC_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCBDDC );
          break;

      case KSP_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCKSP );
          break;

      case PYTHON_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCPYTHON );
          break;

      case AMGX_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCAMGX );
          break;

      case QR_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCQR );
          break;

      case NN_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCNN );
          break;

      case SPAI_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSPAI );
          break;

      case MAT_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCMAT );
          break;

      case HYPRE_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCHYPRE );
          break;

      case PARMS_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCPARMS );
          break;

      case TFS_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCTFS );
          break;

      case GALERKIN_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCGALERKIN );
          break;

      case EXOTIC_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCEXOTIC );
          break;

      case CP_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCCP );
          break;

      case BFBT_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCBFBT );
          break;

      case PFMG_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCPFMG );
          break;

      case SMG_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSMG );
          break;

      case SYSPFMG_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSYSPFMG );
          break;

      case REDISTRIBUTE_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCREDISTRIBUTE );
          break;

      case SVD_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSVD );
          break;

      case CHOWILUVIENNACL_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCCHOWILUVIENNACL );
          break;

      case ROWSCALINGVIENNACL_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCROWSCALINGVIENNACL );
          break;

      case SAVIENNACL_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCSAVIENNACL );
          break;

      case KACZMARZ_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCKACZMARZ );
          break;

      case TELESCOPE_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCTELESCOPE );
          break;

      case PATCH_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCPATCH );
          break;

      case LMVM_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCLMVM );
          break;

      case HMG_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCHMG );
          break;

      case DEFLATION_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCDEFLATION );
          break;

      case HPDDM_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCHPDDM );
          break;

      case H2OPUS_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCH2OPUS );
          break;

      case MPI_PRECOND:
          ierr = PCSetType( M_pc, (char*)PCMPI );
          break;

      case AMG_PRECOND:
          // Historically in Feel++ we had "AMG" to mean generic algebraic MG. 
          // You could do PCGAMG or PCHYPRE
          ierr = PCSetType( M_pc, (char*)PCGAMG );
          break;

      default:
          std::cerr << "ERROR: Unsupported PETSC PreconditionerType: "
                    << this->preconditionerType() << "\n"
                    << "Continuing with PETSc default (none)\n";
          ierr = PCSetType( M_pc, (char*)PCNONE );
          break;
    } // end switch

    CHKERRABORT( this->worldComm().globalComm(), ierr );
} // setPetscPreconditionerType()





//------------------------------------------------------------------
// Explicit instantiations
template class SolverLinearPetsc<double>;

} // Feel

#endif // #ifdef FEELPP_HAS_PETSC

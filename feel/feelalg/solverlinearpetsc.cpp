
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
            std::cout << " " << it  << " " << s->prefix() << " KSP Residual norm " << std::scientific << rnorm << "\n";
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
        preconditioner->init();
        VLOG(2) << "__feel_petsc_preconditioner_setup: init prec\n";
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
            boost::shared_ptr<VectorPetsc<double> > x_vec;
            boost::shared_ptr<VectorPetsc<double> > y_vec;
            if ( preconditioner->worldComm().localSize() > 1 )
            {
                CHECK ( preconditioner->matrix() ) << "matrix is not defined";
                x_vec.reset( new VectorPetscMPI<double>( x, preconditioner->matrix()->mapColPtr() ) );
                y_vec.reset( new VectorPetscMPI<double>( y, preconditioner->matrix()->mapRowPtr() ) );
            }
            else
            {
                x_vec.reset( new VectorPetsc<double>( x ) );
                y_vec.reset( new VectorPetsc<double>( y ) );
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
        bool useInitialGuessNonZero = boption(_name="ksp-use-initial-guess-nonzero", _prefix=this->prefix() );
        ierr = KSPSetInitialGuessNonzero ( M_ksp, (useInitialGuessNonZero)?PETSC_TRUE:PETSC_FALSE );
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
            int nRestartGMRES = ioption(_name="gmres-restart", _prefix=this->prefix() );
            ierr = KSPGMRESSetRestart( M_ksp, nRestartGMRES );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
        }
        else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPFGMRES ) )
        {
            int nRestartFGMRES = ioption(_name="fgmres-restart", _prefix=this->prefix() );
            ierr = KSPGMRESSetRestart( M_ksp, nRestartFGMRES );
            CHKERRABORT( this->worldComm().globalComm(),ierr );
            if ( this->M_preconditioner )
                this->M_preconditioner->setSide( preconditioner_type::RIGHT );
        }
        else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPGCR ) )
        {
            int nRestartGCR = ioption(_name="gcr-restart", _prefix=this->prefix() );
            ierr = KSPGCRSetRestart( M_ksp, nRestartGCR );
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
                                      PETSC_NULL,   // pointer to the array which holds the history
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
#if PETSC_VERSION_LESS_THAN(3,4,0)
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
#if PETSC_VERSION_LESS_THAN(3,4,0)
                KSPSetPreconditionerSide( M_ksp, PC_LEFT );
#else
                KSPSetPCSide( M_ksp, PC_LEFT );
#endif
                break;
            case preconditioner_type::RIGHT:
                VLOG(2) << " . PC is set to right side\n";
#if PETSC_VERSION_LESS_THAN(3,4,0)
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

        if ( Environment::vm(_name="ksp-monitor",_prefix=this->prefix()).template as<bool>() )
        {
            //KSPMonitorSet( M_ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
            KSPMonitorSet( M_ksp,__feel_petsc_monitor,(void*) this,PETSC_NULL );
        }

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
    this->setWorldComm( matrix_in.comm() );
    this->init ();

    MatrixPetsc<T> * matrix   = const_cast<MatrixPetsc<T> *>( dynamic_cast<MatrixPetsc<T> const*>( &matrix_in ) );
    MatrixPetsc<T> * precond  = const_cast<MatrixPetsc<T> *>( dynamic_cast<MatrixPetsc<T> const*>( &precond_in ) );
    VectorPetsc<T> * solution = dynamic_cast<VectorPetsc<T>*>( &solution_in );
    VectorPetsc<T> * rhs      = const_cast<VectorPetsc<T> *>( dynamic_cast<VectorPetsc<T> const*>( &rhs_in ) );

    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
    FEELPP_ASSERT( matrix   != NULL ).error( "non petsc matrix structure" );
    FEELPP_ASSERT( precond  != NULL ).error( "non petsc matrix structure" );
    FEELPP_ASSERT( solution != NULL ).error( "non petsc vector structure" );
    FEELPP_ASSERT( rhs      != NULL ).error( "non petsc vector structure" );

    int ierr=0;
    int its=0;
    PetscReal final_resid=0.;

    // Close the matrices and vectors in case this wasn't already done.
    matrix->close ();
    precond->close ();
    solution->close ();
    rhs->close ();


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

    if ( boption( _prefix=this->prefix(), _name="ksp-view" ) )
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
    double* p;
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
    double* p;
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

        ierr = MatNullSpaceCreate( PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        ierr = KSPSetNullSpace( M_ksp, nullsp );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
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
    ierr = MatNullSpaceCreate( this->worldComm(), PETSC_FALSE , dimNullSpace, petsc_vec.data()/*PETSC_NULL*/, &nullsp );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    //ierr = MatNullSpaceView( nullsp, PETSC_VIEWER_STDOUT_WORLD );
    //CHKERRABORT( this->worldComm().globalComm(),ierr );

    ierr = MatSetNullSpace(A,nullsp);
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    //ierr = MatNullSpaceRemove(nullsp,rhs);
    //CHKERRABORT( this->worldComm().globalComm(),ierr );

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
    ierr = MatNullSpaceCreate( this->worldComm(), PETSC_FALSE , dimNullSpace, petsc_vec.data()/*PETSC_NULL*/, &nullsp );
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
    DVLOG(2) << "[SolverLinearPetsc] solver type:  " << this->solverType() << "\n";

    switch ( this->solverType() )
    {

    case CG:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCG );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case CR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case CGS:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCGS );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case BICG:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPBICG );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case TCQMR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPTCQMR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case TFQMR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPTFQMR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case LSQR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPLSQR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case BICGSTAB:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPBCGS );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case MINRES:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPMINRES );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case GMRES:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPGMRES );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case FGMRES:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPFGMRES );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case RICHARDSON:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPRICHARDSON );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case CHEBYSHEV:
#if PETSC_VERSION_LESS_THAN(3,3,0)
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCHEBYCHEV );
#else
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCHEBYSHEV );
#endif
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case PREONLY :
        ierr = KSPSetType ( M_ksp, ( char* ) KSPPREONLY );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case GCR :
        ierr = KSPSetType ( M_ksp, ( char* ) KSPGCR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Solver: "
                  << this->solverType()               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

}








template <typename T>
void
SolverLinearPetsc<T>::setPetscPreconditionerType()
{

    int ierr = 0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERUMFPACK );

    if ( ierr )
    {
        ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERSUPERLU );

        if ( ierr )
        {
            ierr = PCFactorSetMatSolverPackage( M_pc,MATSOLVERPETSC );
        }
    }

#elif (PETSC_VERSION_MAJOR >= 3)
    ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_UMFPACK );

    if ( ierr )
    {
        ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_SUPERLU );

        if ( ierr )
        {
            ierr = PCFactorSetMatSolverPackage( M_pc,MAT_SOLVER_PETSC );
        }
    }

#endif
    DVLOG(2) << "[SolverLinearPetsc] preconditioner type:  " << this->preconditionerType() << "\n";

    switch ( this->preconditionerType() )
    {
    case IDENTITY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCNONE );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case CHOLESKY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCCHOLESKY );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case ICC_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCICC );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case ILU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCILU );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        if ( this->vm().count( "pc-factor-levels" ) )
        {
            PCFactorSetLevels( M_pc,this->vm()["pc-factor-levels"].template as<int>() );
        }

        else
            PCFactorSetLevels( M_pc,3 );

        if ( this->vm().count( "pc-factor-fill" ) )
        {
            PCFactorSetFill( M_pc,this->vm()["pc-factor-fill"].template as<double>() );
        }

        else
            PCFactorSetFill( M_pc,40 );

        return;

    case LU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCLU );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case ASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCASM );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    case GASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCGASM );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;
#endif

    case JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCJACOBI );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCBJACOBI );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case SOR_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCSOR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case EISENSTAT_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCEISENSTAT );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

#if !((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1) && (PETSC_VERSION_SUBMINOR <= 1))

    case USER_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCMAT );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;
#endif

    case SHELL_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCSHELL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case FIELDSPLIT_PRECOND:
        ierr = PCSetType( M_pc,( char* ) PCFIELDSPLIT );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        ierr = PCFieldSplitSetType( M_pc,PC_COMPOSITE_SCHUR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    case ML_PRECOND:
        ierr = PCSetType( M_pc,( char* ) PCML );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
        return;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << this->preconditionerType()       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

}




//------------------------------------------------------------------
// Explicit instantiations
template class SolverLinearPetsc<double>;

} // Feel

#endif // #ifdef FEELPP_HAS_PETSC

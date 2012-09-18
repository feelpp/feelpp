/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-02

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file solvernonlinearpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-02
 */
#include <feel/feelcore/feel.hpp>



#if defined( FEELPP_HAS_PETSC_H )
#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/solverlinearpetsc.hpp>
#include <feel/feelalg/solvernonlinearpetsc.hpp>
#include <feel/feelalg/functionspetsc.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>


//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
    // Older versions of PETSc do not have the different int typedefs.
    // On 64-bit machines, PetscInt may actually be a long long int.
    // This change occurred in Petsc-2.2.1.
# if (((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 2) && (PETSC_VERSION_SUBMINOR == 0)) || \
      ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)))
    typedef int PetscErrorCode;
    typedef int PetscInt;
#endif

    //-------------------------------------------------------------------
    // this function is called by PETSc at the end of each nonlinear step
    PetscErrorCode
    __feel_petsc_snes_monitor ( SNES snes, PetscInt its, PetscReal fnorm, void *ctx )
    {
        Feel::SolverNonLinearPetsc<double>* solver =
            static_cast<Feel::SolverNonLinearPetsc<double>*> ( ctx );


        //int ierr=0;
        //if (its > 0)
        std::ostringstream ostr;

        ostr << "[SolverNonLinearPetsc] NL step " << its
             << std::scientific
             << ", |residual|_2 = " << fnorm;
        //Feel::Log() << ostr.str() << "\n";
        //std::cout << ostr.str() << "\n";
#if 1
        KSP            ksp;         /* linear solver context */
        SNESGetKSP( snes,&ksp );
        PetscInt lits;
        int ierr = KSPGetIterationNumber ( ksp, &lits );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        PetscReal final_resid;
        // Get the norm of the final residual to return to the user.
        ierr = KSPGetResidualNorm ( ksp, &final_resid );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        ///Feel::Log() << "[SolverNonLinearPetsc] KSP num of it = " << lits << " residual = " << final_resid << "\n";
        //std::cout << "[SolverNonLinearPetsc] KSP num of it = " << lits << " residual = " << final_resid << "\n";
#endif
        KSPConvergedReason reason;
        KSPGetConvergedReason( ksp,&reason );
        if ( reason> 0 )
        {
            if ( solver->showKSPConvergedReason() && solver->worldComm().globalRank() == solver->worldComm().masterRank() && its>0 )
                std::cout<< "  Linear solve converged due to " << Feel::PetscConvertKSPReasonToString(reason)
                         << " iterations " << lits << std::endl;
        }
        else
        {
            if ( solver->showKSPConvergedReason() && solver->worldComm().globalRank() == solver->worldComm().masterRank() && its>0 )
                std::cout<< "  Linear solve did not converge due to " << Feel::PetscConvertKSPReasonToString(reason)
                         << " iterations " << lits << std::endl;
        }


        //return ierr;
        return 0;
    }



    //---------------------------------------------------------------
    // this function is called by PETSc to evaluate the residual at X
    PetscErrorCode
    __feel_petsc_snes_residual ( SNES snes, Vec x, Vec r, void *ctx )
    {
        int ierr=0;

        assert ( x   != NULL );
        assert ( r   != NULL );
        assert ( ctx != NULL );

        Feel::SolverNonLinearPetsc<double>* solver =
            static_cast<Feel::SolverNonLinearPetsc<double>*> ( ctx );

#if !defined(FEELPP_ENABLE_MPI_MODE)
        boost::shared_ptr<Feel::Vector<double> > R( new Feel::VectorPetsc<double>( r ) );
        boost::shared_ptr<Feel::MatrixSparse<double> >  PC;
        boost::shared_ptr<Feel::Vector<double> > X_global( new Feel::VectorPetsc<double>( x ) );
        boost::shared_ptr<Feel::Vector<double> > X_local( new Feel::VectorPetsc<double>( X_global->size() ) );

        X_global->localize ( *X_local );
#else
        boost::shared_ptr<Feel::Vector<double> > R;
        boost::shared_ptr<Feel::Vector<double> > X_global;

        if ( solver->comm().size()>1 )
        {
            R.reset( new Feel::VectorPetscMPI<double>( r,solver->mapRow() ) );
            X_global.reset( new Feel::VectorPetscMPI<double>( x,solver->mapRow() ) );
        }

        else // MPI
        {
            R.reset( new Feel::VectorPetsc<double>( r ) );
            X_global.reset( new Feel::VectorPetsc<double>( x ) );
        }

#endif

        //if (solver->residual != NULL) solver->residual (X_local, R);
        if ( solver->residual != NULL ) solver->residual ( X_global, R );

        //if (solver->matvec   != NULL) solver->matvec   (X_local, R, PC );

        R->close();

        return ierr;
    }



    //---------------------------------------------------------------
    // this function is called by PETSc to evaluate the Jacobian at X
    PetscErrorCode
    __feel_petsc_snes_jacobian ( SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx )
    {
        int ierr=0;
        assert ( ctx != NULL );

        Feel::SolverNonLinearPetsc<double>* solver =
            static_cast<Feel::SolverNonLinearPetsc<double>*> ( ctx );

#if !defined(FEELPP_ENABLE_MPI_MODE)
        boost::shared_ptr<Feel::Vector<double> > R;
        boost::shared_ptr<Feel::MatrixSparse<double> >  PC( new Feel::MatrixPetsc<double>( *pc ) );
        boost::shared_ptr<Feel::MatrixSparse<double> >  Jac( new Feel::MatrixPetsc<double>( *jac ) );
        boost::shared_ptr<Feel::Vector<double> > X_global( new Feel::VectorPetsc<double>( x ) );
        boost::shared_ptr<Feel::Vector<double> > X_local( new Feel::VectorPetsc<double>( X_global->size() ) );

        X_global->localize ( *X_local );
#else // MPI
        boost::shared_ptr<Feel::MatrixSparse<double> > Jac;
        boost::shared_ptr<Feel::Vector<double> > X_global;

        if ( solver->comm().size()>1 )
        {
            Jac.reset( new Feel::MatrixPetscMPI<double>( *jac,solver->mapRow(),solver->mapCol() ) );
            X_global.reset( new Feel::VectorPetscMPI<double>( x,solver->mapRow() ) );
        }

        else
        {
            Jac.reset( new Feel::MatrixPetsc<double>( *jac ) );
            X_global.reset( new Feel::VectorPetsc<double>( x ) );
        }

#endif

        //if (solver->jacobian != NULL) solver->jacobian (X_local, PC );
        if ( solver->jacobian != NULL ) solver->jacobian ( X_global, Jac );

        //if (solver->matvec   != NULL) solver->matvec   (X_local, R, PC );

        //PC=Jac;
        //PC->close();
        //Jac = PC;
        Jac->close();
        //*PC = Jac;
        //PC->printMatlab( "pc.m" );
        //Jac->printMatlab( "jac.m" );
        //!!!*msflag = SAME_NONZERO_PATTERN;
        //*msflag = DIFFERENT_NONZERO_PATTERN;
#if 0
        std::cout << "[__ snes_jacobian] mflag = " << *msflag << "\n";
        std::cout << "[__ snes_jacobian] SAME_NONZERO_PATTERN = " << SAME_NONZERO_PATTERN << "\n";
        std::cout << "[__ snes_jacobian] DIFFERENT_NONZERO_PATTERN = " << DIFFERENT_NONZERO_PATTERN << "\n";
        std::cout << "[__ snes_jacobian] SAME_PRECONDITIONER = " << SAME_PRECONDITIONER << "\n";
#endif

        //std::cout << "\nsolver->reuseJacobian() " << solver->reuseJacobian() << " solver->reusePreconditioner() " << solver->reusePreconditioner() << "\n";
        //solver->setReuse( solver->reuseJacobian(), solver->reusePreconditioner() );
        //solver->setReuse( -1, 2 );
        //*msflag = MatStructure( solver->precMatrixStructure() );
        //*msflag = SAME_NONZERO_PATTERN;
        //*msflag = DIFFERENT_NONZERO_PATTERN;
        return ierr;
    }

    //---------------------------------------------------------------
    // this function is called by PETSc to evaluate the residual at X
    PetscErrorCode
    __feel_petsc_snes_dense_residual ( SNES snes, Vec x, Vec r, void *ctx )
    {
        int ierr=0;

        assert ( x   != NULL );
        assert ( r   != NULL );
        assert ( ctx != NULL );

        Feel::SolverNonLinearPetsc<double>* solver =
            static_cast<Feel::SolverNonLinearPetsc<double>*> ( ctx );

        int size;
        VecGetSize( r,&size );
        double *xa;
        VecGetArray( x, &xa );
        boost::numeric::ublas::vector<double> xx( size );
        std::copy( xa, xa+size, xx.begin() );

        //Feel::Log() << "dense_residual before xx= " << xx << "\n";

        boost::numeric::ublas::vector<double> rr( size );

        if ( solver->dense_residual != NULL ) solver->dense_residual ( xx, rr );

        //Feel::Log() << "dense_residual after update rr= " << rr << "\n";

        for ( int i=0; i<size; i++ )
        {
            VecSetValues( r,1,&i,&rr[i],INSERT_VALUES );
        }

        VecRestoreArray( x, &xa );
        //Feel::Log() << "dense_residual rr= " << rr << "\n";

        return ierr;
    }



    //---------------------------------------------------------------
    // this function is called by PETSc to evaluate the Jacobian at X
    PetscErrorCode
    __feel_petsc_snes_dense_jacobian ( SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx )
    {
        int ierr=0;
        assert ( ctx != NULL );

        Feel::SolverNonLinearPetsc<double>* solver =
            static_cast<Feel::SolverNonLinearPetsc<double>*> ( ctx );

        int size;
        VecGetSize( x,&size );
        double *xa;
        VecGetArray( x, &xa );
        boost::numeric::ublas::vector<double> xx( size );
        std::copy( xa, xa+size, xx.begin() );

        ///Feel::Log() << "dense_jacobian xx= " << xx << "\n";

        int size1;
        int size2;
        MatGetSize( *jac, &size1, &size2 );
        boost::numeric::ublas::matrix<double> jj( size1, size2 );

        if ( solver->dense_jacobian != NULL ) solver->dense_jacobian ( xx, jj );

        //Feel::Log() << "dense_jacobian jj = " << jj << "\n";

        for ( int i = 0; i < size1; ++i )
            for ( int j = 0; j < size2; ++j )
            {
                MatSetValue( *jac, i, j, jj( i, j ), INSERT_VALUES );

            }

        /*
          Assemble matrix
        */
        MatAssemblyBegin( *jac,MAT_FINAL_ASSEMBLY );
        MatAssemblyEnd( *jac,MAT_FINAL_ASSEMBLY );

        VecRestoreArray( x, &xa );

        *msflag = SAME_NONZERO_PATTERN;

        return ierr;
    }

} // end extern "C"
//---------------------------------------------------------------------

namespace Feel
{
// SolverNonLinearPetsc<> methods
template <typename T>
void SolverNonLinearPetsc<T>::clear ()
{
    if ( this->initialized() )
    {
        this->M_is_initialized = false;

        int ierr=0;
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = SNESDestroy( &M_snes );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
#else
        ierr = SNESDestroy( M_snes );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
#endif
    }
}
// SolverNonLinearPetsc<> methods
template <typename T>
void SolverNonLinearPetsc<T>::setReuse ( int jac, int prec )
{

    this->M_reuse_jac=jac;
    this->M_reuse_prec=prec;

    PetscInt cur_jac, cur_prec;
    SNESGetLagJacobian( M_snes, &cur_jac );
    SNESGetLagPreconditioner( M_snes, &cur_prec );
    //std::cout << "[PETSc non linear solver reuse] prev jac=" << cur_jac << " prev prec=" << cur_prec << "\n";

    SNESSetLagJacobian( M_snes, jac );
    SNESSetLagPreconditioner( M_snes, prec );
    SNESGetLagJacobian( M_snes, &cur_jac );
    SNESGetLagPreconditioner( M_snes, &cur_prec );
    //std::cout << "[PETSc non linear solver reuse] prev jac=" << cur_jac << " prev prec=" << cur_prec << "\n";
}
template <typename T>
void SolverNonLinearPetsc<T>::init ()
{
    // Initialize the data structures if not done so already.
    if ( !this->initialized() )
    {
        this->M_is_initialized = true;

        int ierr=0;

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1) && (PETSC_VERSION_SUBMINOR <= 1))

        // At least until Petsc 2.1.1, the SNESCreate had a different calling syntax.
        // The second argument was of type SNESProblemType, and could have a value of
        // either SNES_NONLINEAR_EQUATIONS or SNES_UNCONSTRAINED_MINIMIZATION.
        ierr = SNESCreate( this->worldComm().globalComm(), SNES_NONLINEAR_EQUATIONS, &M_snes );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

#else

        ierr = SNESCreate( this->worldComm().globalComm(),&M_snes );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

#endif
#if PETSC_VERSION_LESS_THAN(2,3,3)
        ierr = SNESSetMonitor ( M_snes, __feel_petsc_snes_monitor,
                                this, PETSC_NULL );
#else
        // API name change in PETSc 2.3.3
        ierr = SNESMonitorSet ( M_snes, __feel_petsc_snes_monitor,
                                this, PETSC_NULL );
#endif
        CHKERRABORT( this->worldComm().globalComm(),ierr );


        ierr = SNESSetFromOptions( M_snes );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

    }


    int ierr=0;

    // if the non linear solver type is define by the user in the code
    switch ( this->getType() )
    {
    case LINE_SEARCH :
    {
        ierr = SNESSetType( M_snes, SNESLS );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }
    break;

    case TRUST_REGION :
    {
        ierr = SNESSetType( M_snes, SNESTR );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }
    break;

    case SELECT_IN_ARGLIST:
        // no-op
        break;
    }


    double __relResTol,__absResTol,__absSolTol;
    int __nbItMax, __nbEvalFuncMax;

    if ( this->getAbsoluteResidualTol()==0 )
    {
        ierr = SNESGetTolerances( M_snes, &__absResTol, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }

    else __absResTol = this->getAbsoluteResidualTol();

    if ( this->getRelativeResidualTol()==0 )
    {
        ierr = SNESGetTolerances( M_snes, PETSC_NULL, &__relResTol, PETSC_NULL, PETSC_NULL, PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }

    else __relResTol = this->getRelativeResidualTol();

    if ( this->getAbsoluteSolutionTol()==0 )
    {
        ierr = SNESGetTolerances( M_snes, PETSC_NULL, PETSC_NULL, &__absSolTol, PETSC_NULL, PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }

    else __absSolTol = this->getAbsoluteSolutionTol();

    if ( this->getNbItMax()==0 )
    {
        ierr = SNESGetTolerances( M_snes, PETSC_NULL, PETSC_NULL, PETSC_NULL, &__nbItMax, PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }

    else __nbItMax = this->getNbItMax();

    ierr = SNESGetTolerances( M_snes, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, &__nbEvalFuncMax );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    ierr = SNESSetTolerances( M_snes,__absResTol,__relResTol,__absSolTol,__nbItMax,__nbEvalFuncMax );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    //KSP ksp;
    ierr = SNESGetKSP ( M_snes, &M_ksp );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    //PC pc;
    ierr = KSPGetPC( M_ksp,&M_pc );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Set user-specified  solver and preconditioner types
    this->setPetscKspSolverType();
    this->setPetscPreconditionerType();
    //this->setPetscConstantNullSpace();
    // sets the software that is used to perform the factorization
    PetscPCFactorSetMatSolverPackage( M_pc,this->matSolverPackageType() );


    if ( this->M_preconditioner )
    {
        //PCSetType( M_pc, PCSHELL );
        PCShellSetContext( M_pc,( void* )this->M_preconditioner.get() );

        //Re-Use the shell functions from petsc_linear_solver
        PCShellSetSetUp( M_pc,__feel_petsc_preconditioner_setup );
        PCShellSetApply( M_pc,__feel_petsc_preconditioner_apply );
    }

    if ( this->showSNESMonitor() )
    {
        ierr = SNESMonitorSet( M_snes,SNESMonitorDefault,PETSC_NULL,PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }

    if ( this->showKSPMonitor() )
    {
        ierr = KSPMonitorSet( M_ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
    }


}



template <typename T>
std::pair< int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearPetsc<T>::solve ( sparse_matrix_ptrtype&  jac_in,  // System Jacobian Matrix
                                 vector_ptrtype& x_in,    // Solution vector
                                 vector_ptrtype& r_in,    // Residual vector
                                 const double,              // Stopping tolerance
                                 const unsigned int )
{
    this->init ();

    int ierr=0;

#if !defined(FEELPP_ENABLE_MPI_MODE)
    MatrixPetsc<T>* jac = dynamic_cast<MatrixPetsc<T>*>( jac_in.get() );
    VectorPetsc<T>* x   = dynamic_cast<VectorPetsc<T>*>( x_in.get() );
    VectorPetsc<T>* r   = dynamic_cast<VectorPetsc<T>*>( r_in.get() );
#else // MPI
    MatrixPetsc<T>* jac;
    VectorPetsc<T>* x;
    VectorPetsc<T>* r;

    if ( this->comm().size()>1 )
    {
        jac = dynamic_cast<MatrixPetscMPI<T>*>( jac_in.get() );
        x = dynamic_cast<VectorPetscMPI<T>*>( x_in.get() );
        r = dynamic_cast<VectorPetscMPI<T>*>( r_in.get() );
        //usefull in __feel_petsc_snes_jacobian and __feel_petsc_snes_residual
        this->setMapRow( jac_in->mapRow() );
        this->setMapCol( jac_in->mapCol() );
    }
    else
    {
        jac = dynamic_cast<MatrixPetsc<T>*>( jac_in.get() );
        x   = dynamic_cast<VectorPetsc<T>*>( x_in.get() );
        r   = dynamic_cast<VectorPetsc<T>*>( r_in.get() );
    }

#endif


    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
    assert( jac != NULL );
    assert( jac->mat() != NULL );
    assert( x   != NULL );
    assert( x->vec()   != NULL );
    assert( r   != NULL );
    assert( r->vec()   != NULL );

    int n_iterations =0;

    ierr = SNESSetFunction ( M_snes, r->vec(), __feel_petsc_snes_residual, this );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    ierr = SNESSetJacobian ( M_snes, jac->mat(), jac->mat(), __feel_petsc_snes_jacobian, this );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    KSPSetOperators( M_ksp, jac->mat(), jac->mat(),
                     MatStructure( ( MatStructure ) this->precMatrixStructure() ) );

    if ( this->preconditionerType() == FIELDSPLIT_PRECOND )
        {
            jac->updatePCFieldSplit( M_pc );
        }

    PreconditionerPetsc<T>::setPetscPreconditionerType( this->preconditionerType(),this->matSolverPackageType(),M_pc, this->worldComm() );

    ierr = SNESSetFromOptions( M_snes );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    //Set the preconditioning matrix
    if ( this->M_preconditioner )
        this->M_preconditioner->setMatrix( jac_in );


    /*
      Set array that saves the function norms.  This array is intended
      when the user wants to save the convergence history for later use
      rather than just to view the function norms via -snes_monitor.
    */
    PetscInt hist_its[50];
    PetscReal history[50];
    ierr = SNESSetConvergenceHistory( M_snes,history,hist_its,50,PETSC_TRUE );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
    // the last one being a pointer to an int to hold the number of iterations required.
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

    ierr = SNESSolve ( M_snes, x->vec(), &n_iterations );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.2.x style
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = SNESSolve ( M_snes, x->vec() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.3.x & newer style
#else

    ierr = SNESSolve ( M_snes, PETSC_NULL, x->vec() );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

#endif

    ierr = SNESGetIterationNumber( M_snes,&n_iterations );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    Log() << "[SolverNonLinearPetsc] number of nonlinear iterations = " << n_iterations << "\n";
    LOG(INFO) << "[SolverNonLinearPetsc] number of nonlinear iterations = " << n_iterations << "\n";

    for ( int i=0; i<n_iterations+1; i++ )
    {
        Log() << "iteration " << i << ": Linear iterations : " << hist_its[i] << " Function norm = " << history[i] << "\n";
        LOG(INFO) << "iteration " << i << ": Linear iterations : " << hist_its[i] << " Function norm = " << history[i] << "\n";
    }

    SNESConvergedReason reason;
    SNESGetConvergedReason( M_snes,&reason );
    Log() << "[solvernonlinearpetsc] convergence reason : " << reason << "\n";
    LOG(INFO) << "[solvernonlinearpetsc] convergence reason : " << reason << "\n";

    CHECK( reason >= 0 ) << "Nonlinear solve did not converge due to " << PetscConvertSNESReasonToString(reason)
                         << " iterations " << n_iterations << std::endl;
    if ( reason<0 )
    {
        Log() << "[solvernonlinearpetsc] not converged: " << reason << "\n";
        LOG(INFO) << "[solvernonlinearpetsc] not converged: " << reason << "\n";
        if (this->showSNESConvergedReason() && this->worldComm().globalRank() == this->worldComm().masterRank() )
            std::cout << "Nonlinear solve did not converge due to " << PetscConvertSNESReasonToString(reason)
                      << " iterations " << n_iterations << std::endl;
    }
    else
    {
        if (this->showSNESConvergedReason() && this->worldComm().globalRank() == this->worldComm().masterRank() )
            std::cout << "Nonlinear solve converged due to " << PetscConvertSNESReasonToString(reason)
                      << " iterations " << n_iterations << std::endl;
    }


#if 0
    this->clear();
#endif
    // return the # of its. and the final residual norm.  Note that
    // n_iterations may be zero for PETSc versions 2.2.x and greater.
    return std::make_pair( reason, 0. );
}


template <typename T>
std::pair<unsigned int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearPetsc<T>::solve ( dense_matrix_type&  jac_in,  // System Jacobian Matrix
                                 dense_vector_type& x_in,    // Solution vector
                                 dense_vector_type& r_in,    // Residual vector
                                 const double,              // Stopping tolerance
                                 const unsigned int )
{
    this->init ();

    int ierr=0;
    Vec petsc_x;
    //VecCreateSeqWithArray(PETSC_COMM_SELF,x_in.size(),x_in.data().data(),&petsc_x);
    VecCreateSeq( this->worldComm().globalComm(),x_in.size(),&petsc_x );

    for ( int i=0; i< ( int )x_in.size(); i++ )
    {
        VecSetValues( petsc_x,1,&i,&x_in[i],INSERT_VALUES );
    }

    Vec petsc_r;
    //VecCreateSeqWithArray(PETSC_COMM_SELF,x_in.size(),x_in.data().data(),&petsc_x);
    VecCreateSeq( this->worldComm().globalComm(),r_in.size(),&petsc_r );

    for ( int i=0; i< ( int )x_in.size(); i++ )
    {
        VecSetValues( petsc_r,1,&i,&r_in[i],INSERT_VALUES );
    }

    ierr = SNESSetFunction ( M_snes, petsc_r, __feel_petsc_snes_dense_residual, this );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    Mat petsc_j;
    MatCreateSeqDense( this->worldComm().globalComm(), jac_in.size1(), jac_in.size2(), 0, &petsc_j );

    for ( int i = 0; i < ( int )jac_in.size1(); ++i )
        for ( int j = 0; j < ( int )jac_in.size2(); ++j )
        {
            MatSetValue( petsc_j, i, j, jac_in( i, j ), INSERT_VALUES );

        }

    /*
      Assemble matrix
    */
    MatAssemblyBegin( petsc_j,MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( petsc_j,MAT_FINAL_ASSEMBLY );


    ierr = SNESSetJacobian ( M_snes, petsc_j, petsc_j, __feel_petsc_snes_dense_jacobian, this );
    CHKERRABORT( this->worldComm().globalComm(),ierr );


    /*
       Set linear solver defaults for this problem. By extracting the
       KSP, KSP, and PC contexts from the SNES context, we can then
       directly call any KSP, KSP, and PC routines to set various options.
    */
    //KSP            ksp;         /* linear solver context */
    //PC             pc;           /* preconditioner context */
    //SNESGetKSP( M_snes,&ksp );
    //KSPGetPC( ksp,&pc );
    //PCSetType(pc,PCNONE);
    PCSetType( M_pc,PCLU );
    KSPSetTolerances( M_ksp,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,20 );

    // Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
    // the last one being a pointer to an int to hold the number of iterations required.

# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
    int n_iterations =0;

    ierr = SNESSolve ( M_snes, petsc_x, &n_iterations );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.2.x style
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = SNESSolve ( M_snes, petsc_x );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // 2.3.x & newer style
#else

    ierr = SNESSolve ( M_snes, PETSC_NULL, petsc_x );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

#endif

    double* a;
    VecGetArray( petsc_x , &a );

    for ( int i = 0; i < ( int )x_in.size(); ++i )
        x_in[i] = a[i];

    VecRestoreArray( petsc_x , &a );

    PETSc::VecDestroy( petsc_x );
    PETSc::VecDestroy( petsc_r );
    PETSc::MatDestroy( petsc_j );

    SNESConvergedReason reason;
    SNESGetConvergedReason( M_snes,&reason );

    //Log() << "[solvernonlinearpetsc] convergence reason : " << reason << "\n";
    if ( reason<0 )
    {
        Debug( 7020 )  << "[solvernonlinearpetsc] not converged (see petscsnes.h for an explanation): " << reason << "\n";
    }

    this->clear();

    // return the # of its. and the final residual norm.  Note that
    // n_iterations may be zero for PETSc versions 2.2.x and greater.
    return std::make_pair( reason, 0. );
}


template <typename T>
void
SolverNonLinearPetsc<T>::setPetscNlSolverType()
{
}

template <typename T>
void
SolverNonLinearPetsc<T>::setPetscKspSolverType()
{
    int ierr = 0;
    Debug( 7010 ) << "[SolverNonLinearPetsc] ksp solver type:  " << this->kspSolverType() << "\n";

    switch ( this->kspSolverType() )
    {

    case CG:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCG );
        CHKERRABORT( this->comm(),ierr );
        return;

    case CR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCR );
        CHKERRABORT( this->comm(),ierr );
        return;

    case CGS:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCGS );
        CHKERRABORT( this->comm(),ierr );
        return;

    case BICG:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPBICG );
        CHKERRABORT( this->comm(),ierr );
        return;

    case TCQMR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPTCQMR );
        CHKERRABORT( this->comm(),ierr );
        return;

    case TFQMR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPTFQMR );
        CHKERRABORT( this->comm(),ierr );
        return;

    case LSQR:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPLSQR );
        CHKERRABORT( this->comm(),ierr );
        return;

    case BICGSTAB:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPBCGS );
        CHKERRABORT( this->comm(),ierr );
        return;

    case MINRES:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPMINRES );
        CHKERRABORT( this->comm(),ierr );
        return;

    case GMRES:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPGMRES );
        CHKERRABORT( this->comm(),ierr );
        return;

    case RICHARDSON:
        ierr = KSPSetType ( M_ksp, ( char* ) KSPRICHARDSON );
        CHKERRABORT( this->comm(),ierr );
        return;

    case CHEBYSHEV:
#if PETSC_VERSION_LESS_THAN(3,3,0)
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCHEBYCHEV );
#else
        ierr = KSPSetType ( M_ksp, ( char* ) KSPCHEBYSHEV );
#endif
        CHKERRABORT( this->comm(),ierr );
        return;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Solver: "
                  << this->kspSolverType()               << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

}

template <typename T>
void
SolverNonLinearPetsc<T>::setPetscPreconditionerType()
{

    int ierr = 0;

    switch ( this->preconditionerType() )
    {
    case IDENTITY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCNONE );
        CHKERRABORT( this->comm(),ierr );
        return;

    case CHOLESKY_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCCHOLESKY );
        CHKERRABORT( this->comm(),ierr );
        return;

    case ICC_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCICC );
        CHKERRABORT( this->comm(),ierr );
        return;

    case ILU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCILU );
        CHKERRABORT( this->comm(),ierr );
        return;

    case LU_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCLU );
        CHKERRABORT( this->comm(),ierr );
        return;

    case ASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCASM );
        CHKERRABORT( this->comm(),ierr );
        return;

    case GASM_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCGASM );
        CHKERRABORT( this->comm(),ierr );
        return;

    case JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCJACOBI );
        CHKERRABORT( this->comm(),ierr );
        return;

    case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCBJACOBI );
        CHKERRABORT( this->comm(),ierr );
        return;

    case SOR_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCSOR );
        CHKERRABORT( this->comm(),ierr );
        return;

    case EISENSTAT_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCEISENSTAT );
        CHKERRABORT( this->comm(),ierr );
        return;

#if !((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1) && (PETSC_VERSION_SUBMINOR <= 1))

    case USER_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCMAT );
        CHKERRABORT( this->comm(),ierr );
        return;
#endif

    case SHELL_PRECOND:
        ierr = PCSetType ( M_pc, ( char* ) PCSHELL );
        CHKERRABORT( this->comm(),ierr );
        return;

    case FIELDSPLIT_PRECOND:
        ierr = PCSetType( M_pc,( char* ) PCFIELDSPLIT );
        CHKERRABORT( this->comm(),ierr );
        return;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << this->preconditionerType()       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

}



//------------------------------------------------------------------
// Explicit instantiations
template class SolverNonLinearPetsc<double>;





} // Feel

#endif // #ifdef FEELPP_HAS_PETSC

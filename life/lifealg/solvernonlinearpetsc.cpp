/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-02

  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file solvernonlinearpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-02
 */
#include <life/lifecore/life.hpp>


#if defined( HAVE_PETSC_H )

#include <life/lifealg/glas.hpp>
#include <life/lifealg/vectorpetsc.hpp>
#include <life/lifealg/matrixpetsc.hpp>
#include <life/lifealg/solvernonlinearpetsc.hpp>


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
  __life_petsc_snes_monitor (SNES snes, PetscInt its, PetscReal fnorm, void *)
  {
    //int ierr=0;

    //if (its > 0)
      std::ostringstream ostr;

      ostr << "[SolverNonLinearPetsc] NL step " << its
           << std::scientific
           << ", |residual|_2 = " << fnorm;
      Life::Log() << ostr.str() << "\n";
#if 1
      KSP            ksp;         /* linear solver context */
      SNESGetKSP( snes,&ksp);
      PetscInt lits;
      int ierr = KSPGetIterationNumber (ksp, &lits);
      CHKERRABORT(Life::Application::COMM_WORLD,ierr);

      PetscReal final_resid;
      // Get the norm of the final residual to return to the user.
      ierr = KSPGetResidualNorm ( ksp, &final_resid);
      CHKERRABORT(Life::Application::COMM_WORLD,ierr);
      Life::Log() << "[SolverNonLinearPetsc] KSP num of it = " << lits << " residual = " << final_resid << "\n";
#endif

    //return ierr;
    return 0;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __life_petsc_snes_residual (SNES snes, Vec x, Vec r, void *ctx)
  {
    int ierr=0;

    assert (x   != NULL);
    assert (r   != NULL);
    assert (ctx != NULL);

    Life::SolverNonLinearPetsc<double>* solver =
      static_cast<Life::SolverNonLinearPetsc<double>*> (ctx);


    boost::shared_ptr<Life::Vector<double> > R( new Life::VectorPetsc<double>(r));
    boost::shared_ptr<Life::MatrixSparse<double> >  PC;
    boost::shared_ptr<Life::Vector<double> > X_global( new Life::VectorPetsc<double>(x));
    boost::shared_ptr<Life::Vector<double> > X_local( new Life::VectorPetsc<double>(X_global->size()) );

    X_global->localize (*X_local);

    //if (solver->residual != NULL) solver->residual (X_local, R);
    if (solver->residual != NULL) solver->residual (X_global, R);
    //if (solver->matvec   != NULL) solver->matvec   (X_local, R, PC );

    R->close();

    return ierr;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __life_petsc_snes_jacobian (SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    int ierr=0;
    assert (ctx != NULL);

    Life::SolverNonLinearPetsc<double>* solver =
      static_cast<Life::SolverNonLinearPetsc<double>*> (ctx);

    boost::shared_ptr<Life::Vector<double> > R;
    boost::shared_ptr<Life::MatrixSparse<double> >  PC(new Life::MatrixPetsc<double>(*pc) );
    boost::shared_ptr<Life::MatrixSparse<double> >  Jac(new Life::MatrixPetsc<double>(*jac) );
    boost::shared_ptr<Life::Vector<double> > X_global( new Life::VectorPetsc<double>(x));
    boost::shared_ptr<Life::Vector<double> > X_local( new Life::VectorPetsc<double>(X_global->size()) );

    X_global->localize (*X_local);

    // if (solver->jacobian != NULL) solver->jacobian (X_local, PC );
    if (solver->jacobian != NULL) solver->jacobian (X_global, Jac );
    //if (solver->matvec   != NULL) solver->matvec   (X_local, R, PC );

    PC->close();
    //Jac = PC;
    Jac->close();

    *msflag = MatStructure( solver->precMatrixStructure() );

    return ierr;
  }

    //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  PetscErrorCode
  __life_petsc_snes_dense_residual (SNES snes, Vec x, Vec r, void *ctx)
  {
    int ierr=0;

    assert (x   != NULL);
    assert (r   != NULL);
    assert (ctx != NULL);

    Life::SolverNonLinearPetsc<double>* solver =
      static_cast<Life::SolverNonLinearPetsc<double>*> (ctx);

    int size;
    VecGetSize(r,&size);
    double *xa; VecGetArray( x, &xa );
    boost::numeric::ublas::vector<double> xx( size );
    std::copy( xa, xa+size, xx.begin() );

    //Life::Log() << "dense_residual before xx= " << xx << "\n";

    boost::numeric::ublas::vector<double> rr( size );

    if (solver->dense_residual != NULL) solver->dense_residual (xx, rr );

    //Life::Log() << "dense_residual after update rr= " << rr << "\n";

    for (int i=0; i<size; i++)
        {
            VecSetValues( r,1,&i,&rr[i],INSERT_VALUES);
        }
    VecRestoreArray( x, &xa );
    //Life::Log() << "dense_residual rr= " << rr << "\n";

    return ierr;
  }



  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  PetscErrorCode
  __life_petsc_snes_dense_jacobian (SNES snes, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    int ierr=0;
    assert (ctx != NULL);

    Life::SolverNonLinearPetsc<double>* solver =
        static_cast<Life::SolverNonLinearPetsc<double>*> (ctx);

    int size;
    VecGetSize(x,&size);
    double *xa; VecGetArray( x, &xa );
    boost::numeric::ublas::vector<double> xx( size );
    std::copy( xa, xa+size, xx.begin() );

    ///Life::Log() << "dense_jacobian xx= " << xx << "\n";

    int size1; int size2;
    MatGetSize( *jac, &size1, &size2 );
    boost::numeric::ublas::matrix<double> jj( size1, size2 );

    if (solver->dense_jacobian != NULL) solver->dense_jacobian (xx, jj );

    //Life::Log() << "dense_jacobian jj = " << jj << "\n";

    for( int i = 0; i < size1; ++i )
        for( int j = 0; j < size2; ++j )
            {
                MatSetValue( *jac, i, j, jj( i, j ), INSERT_VALUES );

            }

    /*
      Assemble matrix
    */
    MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);

    VecRestoreArray( x, &xa );

    *msflag = SAME_NONZERO_PATTERN;

    return ierr;
  }

} // end extern "C"
//---------------------------------------------------------------------

namespace Life
{
// SolverNonLinearPetsc<> methods
template <typename T>
void SolverNonLinearPetsc<T>::clear ()
{
    if (this->initialized())
        {
            this->M_is_initialized = false;

            int ierr=0;

            ierr = SNESDestroy(M_snes);
            CHKERRABORT(Application::COMM_WORLD,ierr);
        }
}

template <typename T>
void SolverNonLinearPetsc<T>::init ()
{
    // Initialize the data structures if not done so already.
    if (!this->initialized())
        {
            this->M_is_initialized = true;

            int ierr=0;

# if ((PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1) && (PETSC_VERSION_SUBMINOR <= 1))

            // At least until Petsc 2.1.1, the SNESCreate had a different calling syntax.
            // The second argument was of type SNESProblemType, and could have a value of
            // either SNES_NONLINEAR_EQUATIONS or SNES_UNCONSTRAINED_MINIMIZATION.
            ierr = SNESCreate(Application::COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &M_snes);
            CHKERRABORT(Application::COMM_WORLD,ierr);

#else

            ierr = SNESCreate(Application::COMM_WORLD,&M_snes);
            CHKERRABORT(Application::COMM_WORLD,ierr);

#endif



            ierr = SNESSetFromOptions(M_snes);
            CHKERRABORT(Application::COMM_WORLD,ierr);
        }
}



template <typename T>
std::pair<unsigned int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearPetsc<T>::solve ( sparse_matrix_ptrtype&  jac_in,  // System Jacobian Matrix
                                 vector_ptrtype& x_in,    // Solution vector
                                 vector_ptrtype& r_in,    // Residual vector
                                 const double,              // Stopping tolerance
                                 const unsigned int)
{
    this->init ();

    int ierr=0;

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR == 3) && (PETSC_VERSION_SUBMINOR >= 3)
    ierr = SNESMonitorSet (M_snes, __life_petsc_snes_monitor, this, PETSC_NULL);
#else
    ierr = SNESSetMonitor (M_snes, __life_petsc_snes_monitor, this, PETSC_NULL);
#endif
    CHKERRABORT(Application::COMM_WORLD,ierr);

    MatrixPetsc<T>* jac = dynamic_cast<MatrixPetsc<T>*>( jac_in.get() );
    VectorPetsc<T>* x   = dynamic_cast<VectorPetsc<T>*>( x_in.get() );
    VectorPetsc<T>* r   = dynamic_cast<VectorPetsc<T>*>( r_in.get() );

    // We cast to pointers so we can be sure that they succeeded
    // by comparing the result against NULL.
    assert(jac != NULL); assert(jac->mat() != NULL);
    assert(x   != NULL); assert(x->vec()   != NULL);
    assert(r   != NULL); assert(r->vec()   != NULL);


    int n_iterations =0;

    ierr = SNESSetFunction (M_snes, r->vec(), __life_petsc_snes_residual, this);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    ierr = SNESSetJacobian (M_snes, jac->mat(), jac->mat(), __life_petsc_snes_jacobian, this);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    /*
      Set array that saves the function norms.  This array is intended
      when the user wants to save the convergence history for later use
      rather than just to view the function norms via -snes_monitor.
    */
    PetscInt hist_its[50];
    PetscReal history[50];
    ierr = SNESSetConvergenceHistory(M_snes,history,hist_its,50,PETSC_TRUE);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
    // the last one being a pointer to an int to hold the number of iterations required.
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

    ierr = SNESSolve (M_snes, x->vec(), &n_iterations);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.2.x style
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = SNESSolve (M_snes, x->vec());
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.3.x & newer style
#else

    ierr = SNESSolve (M_snes, PETSC_NULL, x->vec());
    CHKERRABORT(Application::COMM_WORLD,ierr);

#endif

    ierr = SNESGetIterationNumber(M_snes,&n_iterations);CHKERRABORT(Application::COMM_WORLD,ierr);
    Log() << "[SolverNonLinearPetsc] number of nonlinear iterations = " << n_iterations << "\n";
    for (int i=0; i<n_iterations+1; i++)
        {
            Log() << "iteration " << i << ": Linear iterations : " << hist_its[i] << " Function norm = " << history[i] << "\n";
        }

    SNESConvergedReason reason;
    SNESGetConvergedReason(M_snes,&reason);
    Log() << "[solvernonlinearpetsc] convergence reason : " << reason << "\n";
    if (reason<0) {
        Log() << "[solvernonlinearpetsc] not converged: " << reason << "\n";
    }

    this->clear();

    // return the # of its. and the final residual norm.  Note that
    // n_iterations may be zero for PETSc versions 2.2.x and greater.
    return std::make_pair(n_iterations, 0.);
}


template <typename T>
std::pair<unsigned int, typename SolverNonLinearPetsc<T>::real_type>
SolverNonLinearPetsc<T>::solve ( dense_matrix_type&  jac_in,  // System Jacobian Matrix
                                 dense_vector_type& x_in,    // Solution vector
                                 dense_vector_type& r_in,    // Residual vector
                                 const double,              // Stopping tolerance
                                 const unsigned int)
{
    this->init ();

    int ierr=0;
    int n_iterations =0;

    Vec petsc_x;
    //VecCreateSeqWithArray(PETSC_COMM_SELF,x_in.size(),x_in.data().data(),&petsc_x);
    VecCreateSeq(PETSC_COMM_SELF,x_in.size(),&petsc_x);

    for (int i=0; i< x_in.size(); i++)
        {
            VecSetValues( petsc_x,1,&i,&x_in[i],INSERT_VALUES);
        }

    Vec petsc_r;
    //VecCreateSeqWithArray(PETSC_COMM_SELF,x_in.size(),x_in.data().data(),&petsc_x);
    VecCreateSeq(PETSC_COMM_SELF,r_in.size(),&petsc_r);
    for (int i=0; i<x_in.size(); i++)
        {
            VecSetValues( petsc_r,1,&i,&r_in[i],INSERT_VALUES);
        }

    ierr = SNESSetFunction (M_snes, petsc_r, __life_petsc_snes_dense_residual, this);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    Mat petsc_j;
    MatCreateSeqDense( PETSC_COMM_SELF, jac_in.size1(), jac_in.size2(), 0, &petsc_j );

    for( int i = 0; i < jac_in.size1(); ++i )
        for( int j = 0; j < jac_in.size2(); ++j )
            {
                MatSetValue( petsc_j, i, j, jac_in( i, j ), INSERT_VALUES );

            }

    /*
      Assemble matrix
    */
    MatAssemblyBegin(petsc_j,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_j,MAT_FINAL_ASSEMBLY);


    ierr = SNESSetJacobian (M_snes, petsc_j, petsc_j, __life_petsc_snes_dense_jacobian, this);
    CHKERRABORT(Application::COMM_WORLD,ierr);


    /*
       Set linear solver defaults for this problem. By extracting the
       KSP, KSP, and PC contexts from the SNES context, we can then
       directly call any KSP, KSP, and PC routines to set various options.
    */
    KSP            ksp;         /* linear solver context */
    PC             pc;           /* preconditioner context */
    SNESGetKSP(M_snes,&ksp);
    KSPGetPC(ksp,&pc);
    //PCSetType(pc,PCNONE);
    PCSetType(pc,PCLU);
    KSPSetTolerances(ksp,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,20);

    // Older versions (at least up to 2.1.5) of SNESSolve took 3 arguments,
    // the last one being a pointer to an int to hold the number of iterations required.

# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)

    ierr = SNESSolve (M_snes, petsc_x, &n_iterations);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.2.x style
#elif (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)

    ierr = SNESSolve (M_snes, petsc_x );
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // 2.3.x & newer style
#else

    ierr = SNESSolve (M_snes, PETSC_NULL, petsc_x );
    CHKERRABORT(Application::COMM_WORLD,ierr);

#endif

    double* a;
    VecGetArray( petsc_x , &a );

    for( int i = 0; i < x_in.size(); ++i )
        x_in[i] = a[i];

    VecRestoreArray( petsc_x , &a );

    VecDestroy( petsc_x );
    VecDestroy( petsc_r );
    MatDestroy( petsc_j );

    SNESConvergedReason reason;
    SNESGetConvergedReason(M_snes,&reason);
    //Log() << "[solvernonlinearpetsc] convergence reason : " << reason << "\n";
    if (reason<0) {
        Debug( 7020 )  << "[solvernonlinearpetsc] not converged (see petscsnes.h for an explanation): " << reason << "\n";
    }

    this->clear();

    // return the # of its. and the final residual norm.  Note that
    // n_iterations may be zero for PETSc versions 2.2.x and greater.
    return std::make_pair(n_iterations, 0.);
}




//------------------------------------------------------------------
// Explicit instantiations
template class SolverNonLinearPetsc<double>;





} // Life

#endif // #ifdef HAVE_PETSC

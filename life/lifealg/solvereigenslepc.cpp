/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-04

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file solvereigenslepc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-04
 */
#include <life/lifecore/life.hpp>
#include <life/lifealg/solvereigenslepc.hpp>

namespace Life
{

// must have both slepc and petsc
#if defined(HAVE_SLEPC) && defined(HAVE_PETSC)


/*----------------------- functions ----------------------------------*/
template <typename T>
void
SolverEigenSlepc<T>::clear ()
{
    if (this->initialized())
        {
            this->M_is_initialized = false;

            int ierr=0;

            ierr = EPSDestroy(M_eps);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            // SLEPc default eigenproblem solver
            this->M_eigen_solver_type = ARNOLDI;

        }
}



template <typename T>
void
SolverEigenSlepc<T>::init ()
{
    int ierr=0;

    // Initialize the data structures if not done so already.
    if (!this->initialized())
        {
            this->M_is_initialized = true;

            // Create the eigenproblem solver context
            ierr = EPSCreate (Application::COMM_WORLD, &M_eps);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            ierr = EPSGetIP (M_eps, &M_ip);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            PetscReal eta;
            // Set modified Gram-Schmidt orthogonalization as default
            // and leave other parameters unchanged
#if 0
            EPSOrthogonalizationRefinementType refinement;
            ierr = EPSGetOrthogonalization (M_eps, PETSC_NULL, &refinement, &eta);
            ierr = EPSSetOrthogonalization (M_eps, EPS_MGS_ORTH, refinement, eta);
#else
            IPOrthogonalizationRefinementType refinement;
            ierr = IPGetOrthogonalization (M_ip, PETSC_NULL, &refinement, &eta);
            ierr = IPSetOrthogonalization (M_ip, IP_MGS_ORTH, refinement, eta);
#endif // 0
            CHKERRABORT(Application::COMM_WORLD,ierr);

            // Set user-specified  solver
            setSlepcSolverType();
        }
}



template <typename T>
std::pair<unsigned int, unsigned int>
SolverEigenSlepc<T>::solve (MatrixSparse<T> &matrix_A_in,
                            int nev,                  // number of requested eigenpairs
                            int ncv,                  // number of basis vectors
                            const double tol,         // solver tolerance
                            const unsigned int m_its) // maximum number of iterations
{

    this->init ();

    MatrixPetsc<T>* matrix_A   = dynamic_cast<MatrixPetsc<T>*>(&matrix_A_in);

    int ierr=0;

    // converged eigen pairs and number of iterations
    int nconv=0;
    int its=0;

    // The relative error.
    PetscReal error, re, im;

    // Pointer to vectors of the real parts, imaginary parts.
    PetscScalar kr, ki;

    // Close the matrix and vectors in case this wasn't already done.
    matrix_A->close ();


    // just for debugging, remove this
#if 0
    char mat_file[] = "matA.petsc";
    PetscViewer petsc_viewer;
    ierr = PetscViewerBinaryOpen(Application::COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    ierr = MatView(matrix_A->mat(),petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif

    // Set operators.
    ierr = EPSSetOperators (M_eps, matrix_A->mat(), PETSC_NULL);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    //set the problem type and the position of the spectrum
    setSlepcProblemType();
    setSlepcPositionOfSpectrum();

    // Set eigenvalues to be computed.
    ierr = EPSSetDimensions (M_eps, nev, ncv);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Set the tolerance and maximum iterations.
    ierr = EPSSetTolerances (M_eps, tol, m_its);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Set runtime options, e.g.,
    //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
    // Similar to PETSc, these options will override those specified
    // above as long as EPSSetFromOptions() is called _after_ any
    // other customization routines.
    ierr = EPSSetFromOptions (M_eps);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Solve the eigenproblem.
    ierr = EPSSolve (M_eps);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Get the number of iterations.
    ierr = EPSGetIterationNumber (M_eps, &its);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Get number of converged eigenpairs.
    ierr = EPSGetConverged(M_eps,&nconv);
    CHKERRABORT(Application::COMM_WORLD,ierr);


#if 1 //!defined( NDEBUG )
    // ierr = PetscPrintf(Application::COMM_WORLD,
    //         "\n Number of iterations: %d\n"
    //         " Number of converged eigenpairs: %d\n\n", its, nconv);

    // Display eigenvalues and relative errors.
    ierr = PetscPrintf(Application::COMM_WORLD,
                       "           k           ||Ax-kx||/|kx|\n"
                       "   ----------------- -----------------\n" );
    CHKERRABORT(Application::COMM_WORLD,ierr);

    for(int i=0; i<nconv; i++ )
        {
            ierr = EPSGetEigenpair(M_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            ierr = EPSComputeRelativeError(M_eps, i, &error);
            CHKERRABORT(Application::COMM_WORLD,ierr);

#ifdef USE_COMPLEX_NUMBERS
            re = PetscRealPart(kr);
            im = PetscImaginaryPart(kr);
#else
            re = kr;
            im = ki;
#endif

            if (im != .0)
                {
                    ierr = PetscPrintf(Application::COMM_WORLD," %9f%+9f i %12f\n", re, im, error);
                    CHKERRABORT(Application::COMM_WORLD,ierr);
                }
            else
                {
                    ierr = PetscPrintf(Application::COMM_WORLD,"   %12f       %12f\n", re, error);
                    CHKERRABORT(Application::COMM_WORLD,ierr);
                }
        }

    ierr = PetscPrintf(Application::COMM_WORLD,"\n" );
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif // DEBUG

    // return the number of converged eigenpairs
    // and the number of iterations
    return std::make_pair(nconv, its);

}

template <typename T>
std::pair<unsigned int, unsigned int>
SolverEigenSlepc<T>::solve (MatrixSparse<T> &matrix_A_in,
                            MatrixSparse<T> &matrix_B_in,
                            int nev,                  // number of requested eigenpairs
                            int ncv,                  // number of basis vectors
                            const double tol,         // solver tolerance
                            const unsigned int m_its) // maximum number of iterations
{

    this->init ();

    MatrixPetsc<T>* matrix_A   = dynamic_cast<MatrixPetsc<T>*>(&matrix_A_in);
    MatrixPetsc<T>* matrix_B   = dynamic_cast<MatrixPetsc<T>*>(&matrix_B_in);
    int ierr=0;

    // converged eigen pairs and number of iterations
    int nconv=0;
    int its=0;

    // The relative error.
    PetscReal error, re, im;

    // Pointer to vectors of the real parts, imaginary parts.
    PetscScalar kr, ki;

    // Close the matrix and vectors in case this wasn't already done.
    matrix_A->close ();
    matrix_B->close ();

    // just for debugging, remove this
#if 0
    char mat_file[] = "matA.petsc";
    PetscViewer petsc_viewer;
    ierr = PetscViewerBinaryOpen(Application::COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    ierr = MatView(matrix_A->mat(),petsc_viewer);
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif

    // Set operators.
    ierr = EPSSetOperators (M_eps, matrix_A->mat(), matrix_B->mat());
    //ierr = EPSSetOperators (M_eps, matrix_A->mat(), PETSC_NULL);
    //ierr = EPSSetOperators (M_eps, matrix_B->mat(), PETSC_NULL);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    //set the problem type and the position of the spectrum
    setSlepcProblemType();
    setSlepcPositionOfSpectrum();

    // Set eigenvalues to be computed.
    ierr = EPSSetDimensions (M_eps, nev, ncv);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Set the tolerance and maximum iterations.
    ierr = EPSSetTolerances (M_eps, tol, m_its);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Set runtime options, e.g.,
    //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
    // Similar to PETSc, these options will override those specified
    // above as long as EPSSetFromOptions() is called _after_ any
    // other customization routines.
    ierr = EPSSetFromOptions (M_eps);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Solve the eigenproblem.
    ierr = EPSSolve (M_eps);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Get the number of iterations.
    ierr = EPSGetIterationNumber (M_eps, &its);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    // Get number of converged eigenpairs.
    ierr = EPSGetConverged(M_eps,&nconv);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    /*
      Optional: Get some information from the solver and display it
    */
    EPSGetIterationNumber(M_eps, &its);
    PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);
    int lits;
    EPSGetOperationCounters(M_eps,PETSC_NULL,PETSC_NULL,&lits);
    PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);

    EPSType              type;
    EPSGetType(M_eps,&type);
    PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);

    EPSGetDimensions(M_eps,&nev,PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);

    int maxit;double _tol;
    EPSGetTolerances(M_eps,&_tol,&maxit);
    PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",_tol,maxit);

    PetscScalar shift;
    //STGetShift(_st, &shift );
    PetscPrintf(PETSC_COMM_WORLD," shif=%.4g\n",shift);
#if 1 //!defined( NDEBUG )
    // ierr = PetscPrintf(Application::COMM_WORLD,
    //         "\n Number of iterations: %d\n"
    //         " Number of converged eigenpairs: %d\n\n", its, nconv);

    // Display eigenvalues and relative errors.
    ierr = PetscPrintf(Application::COMM_WORLD,
                       "           k           ||Ax-kx||/|kx|\n"
                       "   ----------------- -----------------\n" );
    CHKERRABORT(Application::COMM_WORLD,ierr);

    for(int i=0; i<nconv; i++ )
        {
            ierr = EPSGetEigenpair(M_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL);
            CHKERRABORT(Application::COMM_WORLD,ierr);

            ierr = EPSComputeRelativeError(M_eps, i, &error);
            CHKERRABORT(Application::COMM_WORLD,ierr);

#ifdef USE_COMPLEX_NUMBERS
            re = PetscRealPart(kr);
            im = PetscImaginaryPart(kr);
#else
            re = kr;
            im = ki;
#endif

            if (im != .0)
                {
                    ierr = PetscPrintf(Application::COMM_WORLD," %9f%+9f i %12e\n", re, im, error);
                    CHKERRABORT(Application::COMM_WORLD,ierr);
                }
            else
                {
                    ierr = PetscPrintf(Application::COMM_WORLD,"   %12f       %12e\n", re, error);
                    CHKERRABORT(Application::COMM_WORLD,ierr);
                }
        }

    ierr = PetscPrintf(Application::COMM_WORLD,"\n" );
    CHKERRABORT(Application::COMM_WORLD,ierr);
#endif // DEBUG

    // return the number of converged eigenpairs
    // and the number of iterations
    return std::make_pair(nconv, its);

}


template <typename T>
void SolverEigenSlepc<T>::setSlepcSolverType()
{
    int ierr = 0;

    switch (this->M_eigen_solver_type)
        {
        case POWER:
            ierr = EPSSetType (M_eps, (char*) EPSPOWER);    CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case SUBSPACE:
            ierr = EPSSetType (M_eps, (char*) EPSSUBSPACE); CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case LAPACK:
            ierr = EPSSetType (M_eps, (char*) EPSLAPACK);   CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case ARNOLDI:
            ierr = EPSSetType (M_eps, (char*) EPSARNOLDI);  CHKERRABORT(Application::COMM_WORLD,ierr); break;

        default:
            std::cerr << "ERROR:  Unsupported SLEPc Eigen Solver: "
                      << this->M_eigen_solver_type         << std::endl
                      << "Continuing with SLEPc defaults" << std::endl;
        }
    EPSType etype;
    ierr = EPSGetType(M_eps,&etype);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    Debug() << "solution method:  " << etype << "\n";
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",etype);
    CHKERRABORT(Application::COMM_WORLD,ierr);
    return;
}

template <typename T>
void
SolverEigenSlepc<T>::setSlepcProblemType()
{
#if 1
    int ierr = 0;

    switch (this->M_eigen_problem_type)
        {
        case NHEP:
            ierr = EPSSetProblemType (M_eps, EPS_NHEP);  CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case GNHEP:
            ierr = EPSSetProblemType (M_eps, EPS_GNHEP); CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case HEP:
            ierr = EPSSetProblemType (M_eps, EPS_HEP);   CHKERRABORT(Application::COMM_WORLD,ierr); break;
        case GHEP:
            ierr = EPSSetProblemType (M_eps, EPS_GHEP);  CHKERRABORT(Application::COMM_WORLD,ierr); break;

        default:
            std::cerr << "ERROR:  Unsupported SLEPc Eigen Problem: "
                      << this->M_eigen_problem_type        << std::endl
                      << "Continuing with SLEPc defaults" << std::endl;
        }
#endif
    return;
}



template <typename T>
void
SolverEigenSlepc<T>:: setSlepcPositionOfSpectrum()
{
    int ierr = 0;

    switch (this->M_position_of_spectrum)
        {
        case LARGEST_MAGNITUDE:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_LARGEST_MAGNITUDE);  CHKERRABORT(Application::COMM_WORLD,ierr); return;
        case SMALLEST_MAGNITUDE:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_SMALLEST_MAGNITUDE); CHKERRABORT(Application::COMM_WORLD,ierr); return;
        case LARGEST_REAL:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_LARGEST_REAL);       CHKERRABORT(Application::COMM_WORLD,ierr); return;
        case SMALLEST_REAL:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_SMALLEST_REAL);      CHKERRABORT(Application::COMM_WORLD,ierr); return;
        case LARGEST_IMAGINARY:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_LARGEST_IMAGINARY);  CHKERRABORT(Application::COMM_WORLD,ierr); return;
        case SMALLEST_IMAGINARY:
            ierr = EPSSetWhichEigenpairs (M_eps, EPS_SMALLEST_IMAGINARY); CHKERRABORT(Application::COMM_WORLD,ierr); return;


        default:
            std::cerr << "ERROR:  Unsupported SLEPc position of spectrum: "
                      << this->M_position_of_spectrum        << std::endl;
            throw std::logic_error( "invalid SLEPc position of spectrum parameter" );
        }
}



template <typename T>
std::pair<typename SolverEigenSlepc<T>::real_type, typename SolverEigenSlepc<T>::real_type>
SolverEigenSlepc<T>::eigenPair(unsigned int i, Vector<T> &solution_in)
{
    int ierr=0;

    PetscReal re, im;

    VectorPetsc<T>* solution = dynamic_cast<VectorPetsc<T>*>(&solution_in);

    // real and imaginary part of the ith eigenvalue.
    PetscScalar kr, ki;

    solution->close();

    ierr = EPSGetEigenpair(M_eps, i, &kr, &ki, solution->vec(), PETSC_NULL);
    CHKERRABORT(Application::COMM_WORLD,ierr);

#ifdef USE_COMPLEX_NUMBERS
    re = PetscRealPart(kr);
    im = PetscImaginaryPart(kr);
#else
    re = kr;
    im = ki;
#endif

    return std::make_pair(re, im);
}


template <typename T>
typename SolverEigenSlepc<T>::real_type
SolverEigenSlepc<T>::relativeError(unsigned int i)
{
    int ierr=0;
    PetscReal error;

    ierr = EPSComputeRelativeError(M_eps, i, &error);
    CHKERRABORT(Application::COMM_WORLD,ierr);

    return error;
}

//------------------------------------------------------------------
// Explicit instantiations
template class SolverEigenSlepc<double>;



#endif // #ifdef HAVE_SLEPC

}

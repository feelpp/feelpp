/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-04

  Copyright (C) 2007, 2009 Université Joseph Fourier (Grenoble I)

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
   \file solvereigenslepc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-04
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelslepc.hpp>
#include <feel/feelalg/solvereigenslepc.hpp>

namespace Feel
{

po::options_description
solvereigenslepc_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";


    //int nev,                  // number of requested eigenpairs
    //int ncv,                  // number of basis vectors
    //const double tol,         // solver tolerance
    //const unsigned int m_its) // maximum number of iterations
    po::options_description _options( "Solver EigenValue Slepc -- " + prefix + " solver options" );
    _options.add_options()
    // solver options
    ( ( _prefix+"slepc-solver-type" ).c_str(), Feel::po::value<std::string>()->default_value( "krylovschur" ), "type of eigenvalue solver" )
    ( ( _prefix+"slepc-nev" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of requested eigenpairs" )
    ( ( _prefix+"slepc-ncv" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "number of basis vectors" )
    ( ( _prefix+"slepc-tol" ).c_str(), Feel::po::value<double>()->default_value( 1e-10 ), "solver tolerance" )
    ( ( _prefix+"slepc-nit" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "maximum number of iterations" );

    return _options;
}

// must have both slepc and petsc
#if defined(FEELPP_HAS_SLEPC) && defined(FEELPP_HAS_PETSC)


/*----------------------- functions ----------------------------------*/
template <typename T>
void
SolverEigenSlepc<T>::clear ()
{
    if ( this->initialized() )
    {
        this->M_is_initialized = false;

        int ierr=0;

        ierr = SLEPc::EPSDestroy( M_eps );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        // SLEPc default eigenproblem solver
        //this->M_eigen_solver_type = ARNOLDI;

    }
}



template <typename T>
void
SolverEigenSlepc<T>::init ()
{
    int ierr=0;

    // Initialize the data structures if not done so already.
    if ( !this->initialized() )
    {
        this->M_is_initialized = true;

        // Create the eigenproblem solver context
        ierr = EPSCreate ( PETSC_COMM_WORLD, &M_eps );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 5)
        ierr = EPSGetBV ( M_eps, &M_ip );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        PetscReal eta;
        // Set modified Gram-Schmidt orthogonalization as default
        // and leave other parameters unchanged
        BVOrthogRefineType refinement;
        ierr = BVGetOrthogonalization ( M_ip, PETSC_NULL, &refinement, &eta );
        ierr = BVSetOrthogonalization ( M_ip, BV_ORTHOG_MGS, refinement, eta );

#else
        ierr = EPSGetIP ( M_eps, &M_ip );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        PetscReal eta;
        // Set modified Gram-Schmidt orthogonalization as default
        // and leave other parameters unchanged
#if 0
        EPSOrthogonalizationRefinementType refinement;
        ierr = EPSGetOrthogonalization ( M_eps, PETSC_NULL, &refinement, &eta );
        ierr = EPSSetOrthogonalization ( M_eps, EPS_MGS_ORTH, refinement, eta );
#else
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 2)
        IPOrthogRefineType refinement;
        ierr = IPGetOrthogonalization ( M_ip, PETSC_NULL, &refinement, &eta );
        ierr = IPSetOrthogonalization ( M_ip, IP_ORTHOG_MGS, refinement, eta );
#elif (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 1)
        IPOrthogonalizationRefinementType refinement;
        ierr = IPGetOrthogonalization ( M_ip, PETSC_NULL, &refinement, &eta );
        ierr = IPSetOrthogonalization ( M_ip, IP_ORTH_MGS, refinement, eta );
#else
        IPOrthogonalizationRefinementType refinement;
        ierr = IPGetOrthogonalization ( M_ip, PETSC_NULL, &refinement, &eta );
        ierr = IPSetOrthogonalization ( M_ip, IP_MGS_ORTH, refinement, eta );
#endif //
#endif // 0

#endif
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        // Set user-specified  solver
        setSlepcSolverType();
    }
}



template <typename T>
typename SolverEigenSlepc<T>::solve_return_type
SolverEigenSlepc<T>::solve ( MatrixSparse<T> &matrix_A_in,
                             int nev,                  // number of requested eigenpairs
                             int ncv,                  // number of basis vectors
                             const double tol,         // solver tolerance
                             const unsigned int m_its ) // maximum number of iterations
{

    this->init ();

    MatrixPetsc<T>* matrix_A   = dynamic_cast<MatrixPetsc<T>*>( &matrix_A_in );

    int ierr=0;

    // converged eigen pairs and number of iterations
    int nconv=0;
    int its=0;

    // Close the matrix and vectors in case this wasn't already done.
    matrix_A->close ();


    // just for debugging, remove this
#if 0
    char mat_file[] = "matA.petsc";
    PetscViewer petsc_viewer;
    ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    ierr = MatView( matrix_A->mat(),petsc_viewer );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif

    // Set operators.
    ierr = EPSSetOperators ( M_eps, matrix_A->mat(), PETSC_NULL );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    //set the problem type and the position of the spectrum
    setSlepcProblemType();
    setSlepcPositionOfSpectrum();
    setSlepcSpectralTransform();

    // Set eigenvalues to be computed.
    ierr = EPSSetDimensions ( M_eps, nev, ncv, 2*ncv );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Set the tolerance and maximum iterations.
    ierr = EPSSetTolerances ( M_eps, this->tolerance(), this->maxIterations() );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Set runtime options, e.g.,
    //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
    // Similar to PETSc, these options will override those specified
    // above as long as EPSSetFromOptions() is called _after_ any
    // other customization routines.
    ierr = EPSSetFromOptions ( M_eps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Solve the eigenproblem.
    ierr = EPSSolve ( M_eps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Get the number of iterations.
    ierr = EPSGetIterationNumber ( M_eps, &its );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Get number of converged eigenpairs.
    ierr = EPSGetConverged( M_eps,&nconv );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );


#if 0 //!defined( NDEBUG )
    // The relative error.
    PetscReal error, re, im;

    // Pointer to vectors of the real parts, imaginary parts.
    PetscScalar kr, ki;


    // ierr = PetscPrintf(PETSC_COMM_WORLD,
    //         "\n Number of iterations: %d\n"
    //         " Number of converged eigenpairs: %d\n\n", its, nconv);

    // Display eigenvalues and relative errors.
    ierr = PetscPrintf( PETSC_COMM_WORLD,
                        "           k           ||Ax-kx||/|kx|\n"
                        "   ----------------- -----------------\n" );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    for ( int i=0; i<nconv; i++ )
    {
        ierr = EPSGetEigenpair( M_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

#if SLEPC_VERSION_LT(3,6,0)
        ierr = EPSComputeRelativeError( M_eps, i, &error );
#else
        ierr = EPSComputeError( M_eps, i, EPS_ERROR_RELATIVE, &error );
#endif
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

#ifdef USE_COMPLEX_NUMBERS
        re = PetscRealPart( kr );
        im = PetscImaginaryPart( kr );
#else
        re = kr;
        im = ki;
#endif

        if ( im != .0 )
        {
            ierr = PetscPrintf( PETSC_COMM_WORLD," %9f%+9f i %12f\n", re, im, error );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }

        else
        {
            ierr = PetscPrintf( PETSC_COMM_WORLD,"   %12e       %12e\n", re, error );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }
    }

    ierr = PetscPrintf( PETSC_COMM_WORLD,"\n" );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif // DEBUG

    // TODO: possible memory leak here
    //VecDestroy( M_mode );
#if PETSC_VERSION_LESS_THAN(3,6,0)
    ierr = MatGetVecs( matrix_A->mat(),PETSC_NULL,&M_mode );
#else
    ierr = MatCreateVecs( matrix_A->mat(),PETSC_NULL,&M_mode );
#endif
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    std::vector<double> ret_error( nconv );

    if ( nconv >= 1 )
    {
#if SLEPC_VERSION_LT(3,6,0)
        ierr = EPSComputeRelativeError( M_eps, nconv, ret_error.data() );
#else
        ierr = EPSComputeError( M_eps, nconv, EPS_ERROR_RELATIVE, ret_error.data() );
#endif
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

    // return the number of converged eigenpairs
    // and the number of iterations
    return boost::make_tuple( nconv, its, ret_error );
}

template <typename T>
typename SolverEigenSlepc<T>::solve_return_type
SolverEigenSlepc<T>::solve ( MatrixSparse<T> &matrix_A_in,
                             MatrixSparse<T> &matrix_B_in,
                             int nev,                  // number of requested eigenpairs
                             int ncv,                  // number of basis vectors
                             const double tol,         // solver tolerance
                             const unsigned int m_its ) // maximum number of iterations
{

    this->init ();


    MatrixPetsc<T>* matrix_A   = dynamic_cast<MatrixPetsc<T>*>( &matrix_A_in );
    MatrixPetsc<T>* matrix_B   = dynamic_cast<MatrixPetsc<T>*>( &matrix_B_in );
    int ierr=0;

    // converged eigen pairs and number of iterations
    int nconv=0;
    int its=0;


    // Close the matrix and vectors in case this wasn't already done.
    matrix_A->close ();
    matrix_B->close ();

    // just for debugging, remove this
#if 0
    char mat_file[] = "matA.petsc";
    PetscViewer petsc_viewer;
    ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, mat_file, PETSC_FILE_CREATE, &petsc_viewer );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    ierr = MatView( matrix_A->mat(),petsc_viewer );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif

    // Set operators.
    ierr = EPSSetOperators ( M_eps, matrix_A->mat(), matrix_B->mat() );
    //ierr = EPSSetOperators (M_eps, matrix_A->mat(), PETSC_NULL);
    //ierr = EPSSetOperators (M_eps, matrix_B->mat(), PETSC_NULL);
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    //set the problem type and the position of the spectrum
    setSlepcProblemType();
    setSlepcPositionOfSpectrum();
    setSlepcSpectralTransform();

    // Set eigenvalues to be computed.
    ierr = EPSSetDimensions ( M_eps, nev, ncv, 2*ncv );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Set the tolerance and maximum iterations.
    ierr = EPSSetTolerances ( M_eps, this->tolerance(), this->maxIterations() );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Set runtime options, e.g.,
    //      -eps_type <type>, -eps_nev <nev>, -eps_ncv <ncv>
    // Similar to PETSc, these options will override those specified
    // above as long as EPSSetFromOptions() is called _after_ any
    // other customization routines.
    ierr = EPSSetFromOptions ( M_eps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Solve the eigenproblem.
    ierr = EPSSolve ( M_eps );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Get the number of iterations.
    ierr = EPSGetIterationNumber ( M_eps, &its );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    // Get number of converged eigenpairs.
    ierr = EPSGetConverged( M_eps,&nconv );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    /*
      Optional: Get some information from the solver and display it
    */
    EPSGetIterationNumber( M_eps, &its );
    //PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR < 5)
    int lits;
    EPSGetOperationCounters( M_eps,PETSC_NULL,PETSC_NULL,&lits );
    //PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);
#endif
#if PETSC_VERSION_LESS_THAN(3,4,0)
    const EPSType              type;
#else
    EPSType              type;
#endif
    EPSGetType( M_eps,&type );
    //PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);

    int ncv1, mdv;
    EPSGetDimensions( M_eps,&nev,&ncv1,&mdv );
    //PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);

    int maxit;
    double _tol;
    EPSGetTolerances( M_eps,&_tol,&maxit );
    //PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",_tol,maxit);

#if 0//!defined( NDEBUG )
    PetscScalar shift;
    //STGetShift(M_st, &shift );
    //PetscPrintf(PETSC_COMM_WORLD," shift=%.4g\n",shift);

    // The relative error.
    PetscReal error, re, im;

    // Pointer to vectors of the real parts, imaginary parts.
    PetscScalar kr, ki;

    // ierr = PetscPrintf(PETSC_COMM_WORLD,
    //         "\n Number of iterations: %d\n"
    //         " Number of converged eigenpairs: %d\n\n", its, nconv);

    // Display eigenvalues and relative errors.
    ierr = PetscPrintf( PETSC_COMM_WORLD,
                        "           k           ||Ax-kx||/|kx|\n"
                        "   ----------------- -----------------\n" );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    for ( int i=0; i<nconv; i++ )
    {
        ierr = EPSGetEigenpair( M_eps, i, &kr, &ki, PETSC_NULL, PETSC_NULL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

#if SLEPC_VERSION_LT(3,6,0)
        ierr = EPSComputeRelativeError( M_eps, i, &error );
#else
        ierr = EPSComputeError( M_eps, i, EPS_ERROR_RELATIVE, &error );
#endif
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

        double norm;
        ierr = EPSComputeResidualNorm( M_eps, i, &norm );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );

#ifdef USE_COMPLEX_NUMBERS
        re = PetscRealPart( kr );
        im = PetscImaginaryPart( kr );
#else
        re = kr;
        im = ki;
#endif

        if ( im != .0 )
        {
            ierr = PetscPrintf( PETSC_COMM_WORLD," %9f%+9f i %12e\n", re, im, error );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }

        else
        {
            ierr = PetscPrintf( PETSC_COMM_WORLD,"   %12e       %12e\n", re, error );
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }

        ierr = PetscPrintf( PETSC_COMM_WORLD, "Residual norm=%12e\n", norm );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
    }

    ierr = PetscPrintf( PETSC_COMM_WORLD,"\n" );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif // DEBUG

    // TODO: possible memory leak here
    //VecDestroy( M_mode );
#if PETSC_VERSION_LESS_THAN(3,6,0)
    ierr = MatGetVecs( matrix_A->mat(),PETSC_NULL,&M_mode );
#else
    ierr = MatCreateVecs( matrix_A->mat(),PETSC_NULL,&M_mode );
#endif
    CHKERRABORT( PETSC_COMM_WORLD,ierr );


    std::vector<double> ret_error( nconv );

    if ( nconv >= 1 )
    {
        for ( int i = 0; i < nconv; ++i )
        {
#if SLEPC_VERSION_LT(3,6,0)
            ierr = EPSComputeRelativeError( M_eps, i, &ret_error[i] );
#else
            ierr = EPSComputeError( M_eps, i, EPS_ERROR_RELATIVE, &ret_error[i] );
#endif
            CHKERRABORT( PETSC_COMM_WORLD,ierr );
        }
    }

    //std::for_each( ret_error.begin(), ret_error.end(), []( double e ) { std::cout << " -- ||A x - lambda B x ||/||x|| = " << e << "\n"; } );
    // return the number of converged eigenpairs
    // and the number of iterations
    return boost::make_tuple( nconv, its, ret_error );
}


template <typename T>
void SolverEigenSlepc<T>::setSlepcSolverType()
{
    int ierr = 0;
#if 1

    switch ( this->M_eigen_solver_type )
    {
    case POWER:
        ierr = EPSSetType ( M_eps, ( char* ) EPSPOWER );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case SUBSPACE:
        ierr = EPSSetType ( M_eps, ( char* ) EPSSUBSPACE );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case LAPACK:
        ierr = EPSSetType ( M_eps, ( char* ) EPSLAPACK );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case ARNOLDI:
        ierr = EPSSetType ( M_eps, ( char* ) EPSARNOLDI );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case KRYLOVSCHUR:
        ierr = EPSSetType ( M_eps, ( char* ) EPSKRYLOVSCHUR );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    case LANCZOS:
        ierr = EPSSetType ( M_eps, ( char* ) EPSLANCZOS );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    default:
        std::cerr << "ERROR:  Unsupported SLEPc Eigen Solver: "
                  << this->M_eigen_solver_type         << std::endl
                  << "Continuing with SLEPc defaults" << std::endl;
    }

#else
    ierr = EPSSetType ( M_eps, ( char* ) EPSKRYLOVSCHUR );
    //ierr = EPSSetType (M_eps, (char*) EPSARPACK);
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif
#if PETSC_VERSION_LESS_THAN(3,4,0)
    const EPSType etype;
#else
    EPSType etype;
#endif
    ierr = EPSGetType( M_eps,&etype );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    VLOG(1) << "solution method:  " << etype << "\n";
    //ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",etype);
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    return;
}

template <typename T>
void
SolverEigenSlepc<T>::setSlepcProblemType()
{
#if 1
    int ierr = 0;

    switch ( this->M_eigen_problem_type )
    {
        // non Hermitian
    case NHEP:
        ierr = EPSSetProblemType ( M_eps, EPS_NHEP );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

        // generalized non-Hermitian
    case GNHEP:
        ierr = EPSSetProblemType ( M_eps, EPS_GNHEP );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

        // Hermitian
    case HEP:
        ierr = EPSSetProblemType ( M_eps, EPS_HEP );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

        // generalized Hermitian
    case GHEP:
        ierr = EPSSetProblemType ( M_eps, EPS_GHEP );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

        // Generalized Non-Hermitian with positive (semi-)deﬁnite B
    case PGNHEP:
        ierr = EPSSetProblemType ( M_eps, EPS_PGNHEP );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        break;

    default:
        std::cerr << "ERROR:  Unsupported SLEPc Eigen Problem: "
                  << this->M_eigen_problem_type        << std::endl
                  << "Continuing with SLEPc defaults" << std::endl;
    }

#endif
    VLOG(1) << "Problem  type:  " <<  this->M_eigen_problem_type  << "\n";
    return;
}



template <typename T>
void
SolverEigenSlepc<T>:: setSlepcPositionOfSpectrum()
{
    int ierr = 0;

    switch ( this->M_position_of_spectrum )
    {
    case LARGEST_MAGNITUDE:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_LARGEST_MAGNITUDE );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case SMALLEST_MAGNITUDE:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_SMALLEST_MAGNITUDE );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case LARGEST_REAL:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_LARGEST_REAL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case SMALLEST_REAL:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_SMALLEST_REAL );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case LARGEST_IMAGINARY:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_LARGEST_IMAGINARY );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;

    case SMALLEST_IMAGINARY:
        ierr = EPSSetWhichEigenpairs ( M_eps, EPS_SMALLEST_IMAGINARY );
        CHKERRABORT( PETSC_COMM_WORLD,ierr );
        return;


    default:
        std::cerr << "ERROR:  Unsupported SLEPc position of spectrum: "
                  << this->M_position_of_spectrum        << std::endl;
        throw std::logic_error( "invalid SLEPc position of spectrum parameter" );
    }

    ierr = EPSSetWhichEigenpairs ( M_eps, EPS_SMALLEST_REAL );
}

template <typename T>
void
SolverEigenSlepc<T>:: setSlepcSpectralTransform()
{
    int ierr = 0;

    ST st;
    ierr = EPSGetST ( M_eps, &st );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    switch ( this->spectralTransform() )
    {

    case SINVERT:
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 1)
        ierr = STSetType( st, STSINVERT );
#else
        ierr = STSetType( st, STSINV );
#endif
        break;

    case FOLD:
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR < 5)
        ierr = STSetType( st, STFOLD );
#else
        CHECK( false ) << "FOLD not supported from slepc 3.5\n";
#endif
        break;

    case CAYLEY:
        ierr = STSetType( st, STCAYLEY );
        break;

    case SHIFT:
    default:
        ierr = STSetType( st, STSHIFT );
        break;
    }

    CHKERRABORT( PETSC_COMM_WORLD,ierr );
}



template <typename T>
typename SolverEigenSlepc<T>::eigenpair_type
SolverEigenSlepc<T>::eigenPair( unsigned int i )
{
    int ierr=0;

    PetscReal re, im;

    // real and imaginary part of the ith eigenvalue.
    PetscScalar kr, ki;

    int s;
    ierr = VecGetSize( M_mode,&s );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    ierr = EPSGetEigenpair( M_eps, i, &kr, &ki, M_mode, PETSC_NULL );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

#ifdef USE_COMPLEX_NUMBERS
    re = PetscRealPart( kr );
    im = PetscImaginaryPart( kr );
#else
    re = kr;
    im = ki;
#endif

    //vector_ptrtype solution( new VectorPetsc<value_type>( s, s ) );
    vector_ptrtype solution;
    if ( this->mapRow().worldComm().globalSize()>1 )
        solution = vector_ptrtype( new VectorPetscMPI<value_type>( M_mode,this->mapRowPtr() ) );
    else
        solution = vector_ptrtype( new VectorPetsc<value_type>( M_mode,this->mapRowPtr() ) );

#if 0
    for ( size_type k = 0; k < solution->map().nLocalDofWithGhost(); ++k )
    {
        std::cout << "sol(k)"<<solution->operator()(k) << std::endl;
    }
#endif
#if 0
    double* a;
    ierr = VecGetArray( M_mode, &a );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    for ( int i = 0; i < s; ++i )
    {
        solution->set( i, a[i] );
    }

    solution->close();

    ierr = VecRestoreArray( M_mode, &a );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
#endif
    return boost::make_tuple( re, im, solution );
}

template <typename T>
typename SolverEigenSlepc<T>::eigenmodes_type
SolverEigenSlepc<T>::eigenModes()
{
    int ierr=0;

    int ncv;
    ierr = EPSGetConverged( M_eps,&ncv );
    CHKERRABORT( PETSC_COMM_WORLD,ierr );
    eigenmodes_type modes;

    for ( int i = 0; i < ncv; ++i )
    {
        eigenpair_type mode = this->eigenPair( i );
        // sort with respect to magnitude
        real_type mod2 = mode.template get<0>()*mode.template get<0>()+mode.template get<1>()*mode.template get<1>();
        modes[mod2] = mode;
    }

    return modes;
}


template <typename T>
typename SolverEigenSlepc<T>::real_type
SolverEigenSlepc<T>::relativeError( unsigned int i )
{
    int ierr=0;
    PetscReal error;

#if SLEPC_VERSION_LT(3,6,0)
    ierr = EPSComputeRelativeError( M_eps, i, &error );
#else
    ierr = EPSComputeError( M_eps, i, EPS_ERROR_RELATIVE, &error );
#endif
    CHKERRABORT( PETSC_COMM_WORLD,ierr );

    return error;
}

//------------------------------------------------------------------
// Explicit instantiations
template class SolverEigenSlepc<double>;



#endif // #ifdef FEELPP_HAS_SLEPC

}

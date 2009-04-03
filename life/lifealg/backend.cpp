/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-12-23

  Copyright (C) 2007-2009 Université Joseph Fourier (Grenoble I)

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
   \file backend.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-12-23
 */
#include <life/lifealg/backend.hpp>
#include <life/lifealg/backendgmm.hpp>
#include <life/lifealg/backendpetsc.hpp>
#include <life/lifealg/backendtrilinos.hpp>

namespace Life
{
template <typename T>
Backend<T>::Backend()
    :
#if defined( HAVE_PETSC_H )
    M_backend    (BACKEND_PETSC),
#else
    M_backend    (BACKEND_GMM),
#endif
    M_nlsolver(),
    M_tolerance( 1e-10 ),
    M_transpose( false ),
    M_maxiter( 1000 )
{
}

template <typename T>
Backend<T>::Backend( Backend const& backend )
    :
    M_backend    ( backend.M_backend ),
    M_nlsolver( backend.M_nlsolver ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_tolerance( backend.M_tolerance ),
    M_transpose( backend.M_transpose ),
    M_maxiter( backend.M_maxiter )
{
}
template <typename T>
Backend<T>::Backend( po::variables_map const& vm, std::string const& prefix )
    :
    M_nlsolver( solvernonlinear_type::build( vm ) ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_tolerance( 1e-10 ),
    M_transpose( false ),
    M_maxiter( 1000 )
{
}
template <typename T>
Backend<T>::~Backend()
{
    //this->clear ();
}
template <typename T>
typename Backend<T>::backend_ptrtype
Backend<T>::build( BackendType bt )
{
    // Build the appropriate solver
    switch ( bt )
        {

        case BACKEND_GMM:
            {
                return backend_ptrtype( new BackendGmm<value_type> );
            }
            break;
#if defined ( HAVE_PETSC_H )
        case BACKEND_PETSC:
            {
                return backend_ptrtype( new BackendPetsc<value_type> );
            }
            break;
#endif

#if defined ( HAVE_TRILINOS_EPETRA )
        case BACKEND_TRILINOS:
            {
                return backend_ptrtype( new BackendTrilinos );
            }
            break;
#endif
        default:
            std::cerr << "ERROR:  Unrecognized backend type package: "
                      << bt
                      << std::endl;
            throw std::invalid_argument( "invalid backend type" );
        }

    return backend_ptrtype();
}
template <typename T>
typename Backend<T>::backend_ptrtype
Backend<T>::build( po::variables_map const& vm, std::string const& prefix )
{
    Log() << "[Backend] backend " << vm["backend"].template as<std::string>() << "\n";
    BackendType bt;
    if ( vm["backend"].template as<std::string>() == "gmm" )
        bt = BACKEND_GMM;
    else if ( vm["backend"].template as<std::string>() == "petsc" )
        bt = BACKEND_PETSC;
    else if ( vm["backend"].template as<std::string>() == "trilinos" )
        bt = BACKEND_TRILINOS;
    else
        {

#if defined( HAVE_PETSC_H )

            Log() << "[Backend] use fallback backend petsc\n";
            bt = BACKEND_PETSC;
#else
            Log() << "[Backend] backend " << vm["backend"].template as<std::string>() << " not available\n";
            Log() << "[Backend] use fallback backend gmm\n";
            bt = BACKEND_GMM;
#endif
        }

    // Build the appropriate solver
    switch ( bt )
        {
#if defined ( HAVE_PETSC_H )
        case BACKEND_PETSC:
            {
                Log() << "[Backend] Instantiate a Petsc backend\n";
                return backend_ptrtype( new BackendPetsc<value_type>( vm, prefix ) );
            }
            break;
#endif
#if defined ( HAVE_TRILINOS_EPETRA )
        case BACKEND_TRILINOS:
            {
#if defined ( HAVE_TRILINOS_EPETRA )
                return backend_ptrtype( new BackendTrilinos( vm, prefix ) );
#else
                return backend_ptrtype();
#endif
            }
            break;
#endif
        case BACKEND_GMM:
        default:
            {
                return backend_ptrtype( new BackendGmm<value_type>( vm, prefix ) );
            }
            break;
        }

    // should never happen
    return backend_ptrtype();
}
template <typename T>
typename Backend<T>::solve_return_type
Backend<T>::solve( sparse_matrix_ptrtype const& A,
                   sparse_matrix_ptrtype const& P,
                   vector_ptrtype& x,
                   vector_ptrtype const& b,
                   bool reusePC )
{
    M_reusePC = reusePC;
    if ( !M_reusePC ) {
        reset();
    }
    start();

    this->setPrecMatrixStructure( SAME_PRECONDITIONER );
    boost::tie( M_converged, M_iteration, M_residual ) = this->solve( A, P, x, b );

    stop();
    M_reuseFailed = M_reusedPC && (!M_converged);
    if ( M_reuseFailed )
        {
            reset();
            start();
            this->setPrecMatrixStructure( SAME_NONZERO_PATTERN );
            boost::tie( M_converged, M_iteration, M_residual ) = this->solve( A, P, x, b );
            stop();
        }

    return boost::make_tuple( M_converged, M_iteration, M_residual );
}
template <typename T>
typename Backend<T>::nl_solve_return_type
Backend<T>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its,
                     bool reusePC )
{
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );
    M_nlsolver->solve( A, x, b, tol, its );
    return boost::make_tuple( true, its, tol );

}
template <typename T>
typename Backend<T>::nl_solve_return_type
Backend<T>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its )
{
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );
    M_nlsolver->solve( A, x, b, tol, its );
    return boost::make_tuple( true, its, tol );
}
template <typename T>
typename Backend<T>::real_type
Backend<T>::dot( vector_type const& x, vector_type const& y ) const
{
    real_type localres = 0;
    for( size_type i = 0; i < x.localSize(); ++i )
        {
            localres += x(i)*y(i);
        }
    real_type globalres=localres;
    mpi::all_reduce( Application::comm(), localres, globalres, std::plus<real_type>() );
    return globalres;
}

template <typename T>
void
Backend<T>::start()
{
    M_timer.restart();
}

template <typename T>
void
Backend<T>::stop()
{
    double solveTime = M_timer.elapsed();
    double solveIter = M_iteration + 0.01;
    M_reusedPC = M_reusePC;
    ++M_nUsePC;
    if ( M_nUsePC == 1 )
        {
            M_reusePC = true;
            M_firstSolveTime = solveTime;
            if ( !M_reuseFailed )
                M_maxiter = std::min( M_maxiter, (size_type)(1.5*solveIter + 10.5));
        }
    else
        {
            double nextSolveIter;
            if ( M_nUsePC == 2 )
                {
                    M_totalSolveIter = solveIter*(1.0+M_firstSolveTime/solveTime);
                    nextSolveIter = solveIter;
                }
            else
                {
                    M_totalSolveIter += solveIter;
                    //                 if ( solveIter > M_lastSolveIter )
                    //                     nextSolveIter = 2*solveIter - M_lastSolveIter;
                    //                 else
                    //                     nextSolveIter = solveIter * solveIter / M_lastSolveIter;
                    nextSolveIter = solveIter;
                }
            M_reusePC = ( M_totalSolveIter > M_nUsePC * nextSolveIter );
            M_lastSolveIter = solveIter;
            if ( M_reusePC )
                {
                    M_maxiter = std::min( M_maxiter, (size_type)( M_totalSolveIter/M_nUsePC + 0.5 ) );
                }
        }
}

template <typename T>
void
Backend<T>::reset()
{
    M_reusePC = false;
    M_totalSolveIter = 0.0;
    M_nUsePC = 0;
    //M_backend->set_maxiter( M_maxiter );

}
/*
 * Explicit instantiations
 */
template class Backend<double>;

/**
 * \return the command lines options of the petsc backend
 */
po::options_description backend_options()
{
    po::options_description _options( "Linear and NonLinear Solvers Backend options");
    _options.add_options()
        // solver options
#if defined( HAVE_PETSC_H )
        ("backend", Life::po::value<std::string>()->default_value( "petsc" ), "backend type: gmm, petsc, trilinos")
#else
        ("backend", Life::po::value<std::string>()->default_value( "gmm" ), "backend type: gmm, petsc, trilinos")
#endif
        ;
    return _options;
}

}

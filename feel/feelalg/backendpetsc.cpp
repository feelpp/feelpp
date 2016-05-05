/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-30

  Copyright (C) 2007, 2009 Universit� Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
   \file backendpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-30
 */
#define FEELPP_BACKEND_PETSC_NOEXTERN 1

#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feelcore/feel.hpp>

#if defined( FEELPP_HAS_PETSC_H )

extern "C" {
#include <petscmat.h>
}

namespace Feel
{

PetscErrorCode feel_petsc_post_solve( KSP ksp, Vec x, Vec y, void* ctx )
{
    BackendPetsc<double>* b = static_cast<BackendPetsc<double>*>( ctx );
    LOG( INFO ) << "call feel_petsc_post_solve";

    if ( b->postSolve() )
    {
        vector_ptrtype vx( vec( x, b->dataMap() ) );
        vector_ptrtype vy( vec( y, b->dataMap() ) );
        b->postSolve( vx, vy );
    }
    else
    {
        LOG( WARNING ) << "call feel_petsc_post_solve without post solve function, we duplicate the input into the output";
        auto e = VecDuplicate( x, &y );
        CHKERRABORT( b->comm().globalComm(), e );
        e = VecCopy( x, y );
        CHKERRABORT( b->comm().globalComm(), e );
    }
    LOG( INFO ) << "call feel_petsc_post_solve done";
    return 0;
}
PetscErrorCode feel_petsc_pre_solve( KSP ksp, Vec x, Vec y, void* ctx )
{
    BackendPetsc<double>* b = static_cast<BackendPetsc<double>*>( ctx );

    if ( b->preSolve() )
    {
        vector_ptrtype vx( vec( x, b->dataMap() ) );
        vector_ptrtype vy( vec( y, b->dataMap() ) );
        b->preSolve( vx, vy );
    }
    else
    {
        LOG( WARNING ) << "call feel_petsc_pre_solve without pre solve function, we duplicate the input into the output";
        auto e = VecDuplicate( x, &y );
        CHKERRABORT( b->comm().globalComm(), e );
        e = VecCopy( x, y );
        CHKERRABORT( b->comm().globalComm(), e );
    }
    LOG( INFO ) << "call feel_petsc_pre_solve";
    return 0;
}

template <typename T>
BackendPetsc<T>::~BackendPetsc()
{
    this->clear();
}
template <typename T>
void BackendPetsc<T>::clear()
{
    LOG( INFO ) << "Deleting linear solver petsc";
    M_solver_petsc.clear();
    //LOG(INFO) << "Deleting non linear solver petsc";
    //M_nl_solver_petsc.clear();
    LOG( INFO ) << "Deleting backend petsc";

    super::clear();
}

template <typename T>
typename BackendPetsc<T>::vector_ptrtype
BackendPetsc<T>::toBackendVectorPtr( vector_type const& v )
{
    typedef VectorUblas<T> vector_ublas_type;
    typedef typename vector_ublas_type::range::type vector_ublas_range_type;
    typedef typename vector_ublas_type::slice::type vector_ublas_slice_type;

    vector_ublas_type* vecUblas = const_cast<vector_ublas_type*>( dynamic_cast<vector_ublas_type const*>( &v ) );
    if ( vecUblas )
    {
        //std::cout << "Convert Ublas vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblas ) );
        return _newvec;
    }
    //vector_ublas_range_type * vecUblasRange = dynamic_cast<vector_ublas_range_type*>( &v );
    vector_ublas_range_type* vecUblasRange = const_cast<vector_ublas_range_type*>( dynamic_cast<vector_ublas_range_type const*>( &v ) );
    if ( vecUblasRange )
    {
        //std::cout << "Convert Ublas range vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblasRange ) );
        return _newvec;
    }
#if 0
    vector_ublas_slice_type * vecUblasSlice = dynamic_cast<vector_ublas_slice_type*>( &v );
    if ( vecUblasSlice )
    {
        //std::cout << "Convert Ublas slice vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblasSlice ) );
        return _newvec;
    }
#endif

    //std::cout << "DefaultConvert vector\n";
    return super::toBackendVectorPtr( v );
}

template <typename T>
typename BackendPetsc<T>::vector_ptrtype
BackendPetsc<T>::toBackendVectorPtr( vector_ptrtype const& v )
{
    if ( this->comm().globalSize() > 1 )
    {
        petscMPI_vector_ptrtype vecPetsc = boost::dynamic_pointer_cast<petscMPI_vector_type>( v );
        //if ( vecPetsc ) std::cout << "Convert PetscMPI vector\n";
        if ( vecPetsc )
            return v;
    }
    else
    {
        petsc_vector_ptrtype vecPetsc = boost::dynamic_pointer_cast<petsc_vector_type>( v );
        //if ( vecPetsc ) std::cout << "Convert Petsc vector\n";
        if ( vecPetsc )
            return v;
    }

    typedef VectorUblas<T> vector_ublas_type;
    typedef typename vector_ublas_type::range::type vector_ublas_range_type;
    typedef typename vector_ublas_type::slice::type vector_ublas_slice_type;
    boost::shared_ptr<vector_ublas_type> vecUblas = boost::dynamic_pointer_cast<vector_ublas_type>( v );
    if ( vecUblas )
    {
        //std::cout << "Convert Ublas vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblas ) );
        return _newvec;
    }
    boost::shared_ptr<vector_ublas_range_type> vecUblasRange = boost::dynamic_pointer_cast<vector_ublas_range_type>( v );
    if ( vecUblasRange )
    {
        //std::cout << "Convert Ublas range vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblasRange ) );
        return _newvec;
    }
#if 0
    boost::shared_ptr<vector_ublas_slice_type> vecUblasSlice = boost::dynamic_pointer_cast< vector_ublas_slice_type >( v );
    if ( vecUblasSlice )
    {
        //std::cout << "Convert Ublas slice vector\n";
        vector_ptrtype _newvec( toPETScPtr( *vecUblasSlice ) );
        return _newvec;
    }
#endif

    //std::cout << "DefaultConvert vector\n";
    return super::toBackendVectorPtr( v );
}

template <typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_ptrtype const& A,
                        sparse_matrix_ptrtype const& B,
                        vector_ptrtype& x,
                        vector_ptrtype const& b )
{
    M_solver_petsc.setPrefix( this->prefix() );
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if ( !M_solver_petsc.initialized() )
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.attachNullSpace( this->M_nullSpace );
    M_solver_petsc.attachNearNullSpace( this->M_nearNullSpace );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance = this->rTolerance(),
                                  _atolerance = this->aTolerance(),
                                  _dtolerance = this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    if ( this->preSolve() )
    {
        int e;
        e = KSPSetPreSolve( M_solver_petsc.ksp(), feel_petsc_pre_solve, this );
        CHKERRABORT( this->comm().globalComm(), e );
    }

    if ( this->postSolve() )
    {
        int e;
        e = KSPSetPostSolve( M_solver_petsc.ksp(), feel_petsc_post_solve, this );
        CHKERRABORT( this->comm().globalComm(), e );
    }

    auto res = M_solver_petsc.solve( *A, *B, *x, *b, this->rTolerance(), this->maxIterations(), this->transpose() );
    DVLOG( 2 ) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    DVLOG( 2 ) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.template get<0>() )
        LOG( ERROR ) << "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve

template <typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_type const& A,
                        vector_type& x,
                        vector_type const& b )
{
    M_solver_petsc.setPrefix( this->prefix() );
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if ( !M_solver_petsc.initialized() )
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.attachNullSpace( this->M_nullSpace );
    M_solver_petsc.attachNearNullSpace( this->M_nearNullSpace );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance = this->rTolerance(),
                                  _atolerance = this->aTolerance(),
                                  _dtolerance = this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    if ( this->preSolve() )
    {
        int e;
        e = KSPSetPreSolve( M_solver_petsc.ksp(), feel_petsc_pre_solve, this );
        CHKERRABORT( this->comm().globalComm(), e );
    }

    if ( this->postSolve() )
    {
        int e;
        e = KSPSetPostSolve( M_solver_petsc.ksp(), feel_petsc_post_solve, this );
        CHKERRABORT( this->comm().globalComm(), e );
    }
    auto res = M_solver_petsc.solve( A, x, b, this->rTolerance(), this->maxIterations(), false );
    DVLOG( 2 ) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    DVLOG( 2 ) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.template get<0>() )
        LOG( ERROR ) << "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve

template <typename T>
int BackendPetsc<T>::PtAP( sparse_matrix_ptrtype const& A_,
                           sparse_matrix_ptrtype const& P_,
                           sparse_matrix_ptrtype& C_ ) const
{
    A_->PtAP( *P_, *C_ );
    return 1;
}

template <typename T>
int BackendPetsc<T>::PAPt( sparse_matrix_ptrtype const& A_,
                           sparse_matrix_ptrtype const& P_,
                           sparse_matrix_ptrtype& C_ ) const
{
    A_->PAPt( *P_, *C_ );
    return 1;
}

template <typename T>
void BackendPetsc<T>::prod( sparse_matrix_type const& A,
                            vector_type const& x,
                            vector_type& b, bool transpose ) const
{
    A.multVector( x, b, transpose );
}

template <typename T>
int BackendPetsc<T>::diag( sparse_matrix_type const& A_,
                           vector_type& d_ ) const
{
    A_.diagonal( d_ );
    return 1;
}

template <typename T>
int BackendPetsc<T>::diag( vector_type const& d_,
                           sparse_matrix_type& A_ ) const
{
    A_.setDiagonal( d_ );
    return 1;
}

/**
 * \return the command lines options of the petsc backend
 */
po::options_description backendpetsc_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "BackendPetsc " + prefix + " solver options" );
    _options.add_options()
        // solver options
        ( ( _prefix + "petsc-solver-type" ).c_str(), Feel::po::value<std::string>()->default_value( "umfpack" ), "umfpack, superlu, cg, bicgstab, gmres" )

        // preconditioner options
        ( ( _prefix + "petsc-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "ilut" ), "ilut, ilutp, diag, id" )( ( _prefix + "petsc-threshold" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "threshold value for preconditioners" )( ( _prefix + "petsc-fillin" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "fill-in level value for preconditioners" )

        // solver control options
        ( ( _prefix + "petsc-restart" ).c_str(), Feel::po::value<int>()->default_value( 20 ), "number of iterations before solver restarts (gmres)" )( ( _prefix + "petsc-verbose" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "(=0,1,2) print solver iterations" )( ( _prefix + "petsc-maxiter" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "set maximum number of iterations" )( ( _prefix + "petsc-tolerance" ).c_str(), Feel::po::value<double>()->default_value( 2e-10 ), "set solver tolerance" );
    return _options;
}

/*
 * Explicit instantiations
 */
template class BackendPetsc<double>;
//template class BackendPetsc<std::complex<double>>;

} // Feel
#endif /* FEELPP_HAS_PETSC_H */

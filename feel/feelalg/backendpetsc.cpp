/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-30

  Copyright (C) 2007, 2009 Universit� Joseph Fourier (Grenoble I)

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

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backendpetsc.hpp>

#if defined( FEELPP_HAS_PETSC_H )

extern "C"
{
#include <petscmat.h>
}

namespace Feel
{


template<typename T>
BackendPetsc<T>::~BackendPetsc()
{
    this->clear();
}
template<typename T>
void
BackendPetsc<T>::clear()
{
    LOG(INFO) << "Deleting linear solver petsc";
    M_solver_petsc.clear();
    //LOG(INFO) << "Deleting non linear solver petsc";
    //M_nl_solver_petsc.clear();
    LOG(INFO) << "Deleting backend petsc";

    super::clear();

}

template<typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_ptrtype const& A,
                        sparse_matrix_ptrtype const& B,
                        vector_ptrtype& x,
                        vector_ptrtype const& b )
{
    M_solver_petsc.setPrefix( this->prefix() );
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if (!M_solver_petsc.initialized())
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.attachNullSpace( this->M_nullSpace );
    M_solver_petsc.attachNearNullSpace( this->M_nearNullSpace );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance=this->rTolerance(),
                                  _atolerance=this->aTolerance(),
                                  _dtolerance=this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    auto res = M_solver_petsc.solve( *A, *B, *x, *b, this->rTolerance(), this->maxIterations(), this->transpose() );
    DVLOG(2) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    DVLOG(2) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.template get<0>() )
        LOG(ERROR) << "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve


template<typename T>
typename BackendPetsc<T>::solve_return_type
BackendPetsc<T>::solve( sparse_matrix_type const& A,
                        vector_type& x,
                        vector_type const& b )
{
    M_solver_petsc.setPrefix( this->prefix() );
    M_solver_petsc.setPreconditionerType( this->pcEnumType() );
    M_solver_petsc.setSolverType( this->kspEnumType() );
    if (!M_solver_petsc.initialized())
        M_solver_petsc.attachPreconditioner( this->M_preconditioner );
    M_solver_petsc.setConstantNullSpace( this->hasConstantNullSpace() );
    M_solver_petsc.attachNullSpace( this->M_nullSpace );
    M_solver_petsc.attachNearNullSpace( this->M_nearNullSpace );
    M_solver_petsc.setFieldSplitType( this->fieldSplitEnumType() );
    M_solver_petsc.setTolerances( _rtolerance=this->rTolerance(),
                                  _atolerance=this->aTolerance(),
                                  _dtolerance=this->dTolerance(),
                                  _maxit = this->maxIterations() );
    M_solver_petsc.setPrecMatrixStructure( this->precMatrixStructure() );
    M_solver_petsc.setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_solver_petsc.setShowKSPMonitor( this->showKSPMonitor() );
    M_solver_petsc.setShowKSPConvergedReason( this->showKSPConvergedReason() );

    auto res = M_solver_petsc.solve( A, x, b, this->rTolerance(), this->maxIterations(), false );
    DVLOG(2) << "[BackendPetsc::solve] number of iterations : " << res.template get<1>() << "\n";
    DVLOG(2) << "[BackendPetsc::solve]             residual : " << res.template get<2>() << "\n";

    if ( !res.template get<0>() )
        LOG(ERROR) << "Backend " << this->prefix() << " : linear solver failed to converge" << std::endl;

    return res;
} // BackendPetsc::solve

template <typename T>
int
BackendPetsc<T>::PtAP( sparse_matrix_ptrtype const& A_,
                       sparse_matrix_ptrtype const& P_,
                       sparse_matrix_ptrtype & C_ ) const
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 3)
    MatrixPetsc<T> const* A = dynamic_cast<MatrixPetsc<T> const*> ( A_.get() );
    MatrixPetsc<T> const* P = dynamic_cast<MatrixPetsc<T> const*> ( P_.get() );
    MatrixPetsc<T>* C = dynamic_cast<MatrixPetsc<T>*> ( C_.get() );
    
    return MatPtAP( A->mat(), P->mat(), MAT_INITIAL_MATRIX, 1.0, &C->mat() );
#else
    return -1;
#endif
}


template <typename T>
int
BackendPetsc<T>::PAPt( sparse_matrix_ptrtype const& A_,
                       sparse_matrix_ptrtype const& P_,
                       sparse_matrix_ptrtype & C_ ) const
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 3)
    MatrixPetsc<T> const* A = dynamic_cast<MatrixPetsc<T> const*> ( A_.get() );
    MatrixPetsc<T> const* P = dynamic_cast<MatrixPetsc<T> const*> ( P_.get() );
    MatrixPetsc<T>* C = dynamic_cast<MatrixPetsc<T>*> ( C_.get() );
    
    return MatRARt( A->mat(), P->mat(), MAT_INITIAL_MATRIX, 1.0, &C->mat() );
#else
    return -1;
#endif
}

template <typename T>
int
BackendPetsc<T>::diag( sparse_matrix_type const& A_,
                       vector_type& d_ ) const
{
    MatrixPetsc<T> const& A = dynamic_cast<MatrixPetsc<T>const &> ( A_ );
    VectorPetsc<T> & d = dynamic_cast<VectorPetsc<T>&> ( d_);
    return MatGetDiagonal(A.mat(),d.vec());
}

template <typename T>
int
BackendPetsc<T>::diag( vector_type const & d_,
                       sparse_matrix_type& A_ ) const
{
    MatrixPetsc<T> & A = dynamic_cast<MatrixPetsc<T> &> ( A_ );
    VectorPetsc<T> const& d = dynamic_cast<VectorPetsc<T> const&> ( d_ );
    return MatDiagonalSet(A.mat(),d.vec(),ADD_VALUES);
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
    ( ( _prefix+"petsc-solver-type" ).c_str(), Feel::po::value<std::string>()->default_value( "umfpack" ), "umfpack, superlu, cg, bicgstab, gmres" )

    // preconditioner options
    ( ( _prefix+"petsc-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "ilut" ), "ilut, ilutp, diag, id" )
    ( ( _prefix+"petsc-threshold" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "threshold value for preconditioners" )
    ( ( _prefix+"petsc-fillin" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "fill-in level value for preconditioners" )

    // solver control options
    ( ( _prefix+"petsc-restart" ).c_str(), Feel::po::value<int>()->default_value( 20 ), "number of iterations before solver restarts (gmres)" )
    ( ( _prefix+"petsc-verbose" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "(=0,1,2) print solver iterations" )
    ( ( _prefix+"petsc-maxiter" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "set maximum number of iterations" )
    ( ( _prefix+"petsc-tolerance" ).c_str(), Feel::po::value<double>()->default_value( 2e-10 ), "set solver tolerance" )
    ;
    return _options;
}




/*
 * Explicit instantiations
 */
template class BackendPetsc<double>;
//template class BackendPetsc<std::complex<double>>;



} // Feel
#endif /* FEELPP_HAS_PETSC_H */

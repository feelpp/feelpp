/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-12-23

  Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-present Feel++ Consortium

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
   \file backend.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-12-23
 */
#define FEELPP_BACKEND_NOEXTERN 1

#include <boost/serialization/complex.hpp>

#include <feel/feelalg/backend.hpp>
// PETSc defines MatType which is used as a typename by eigen3 and it conflicts
// undef MatType here to ensure compilation since it is not needed to compile this file
#undef MatType
#include <feel/feelalg/backendeigen.hpp>
#if FEELPP_HAS_PETSC
#include <feel/feelalg/backendpetsc.hpp>
#endif
#include <feel/feelalg/backendtrilinos.hpp>


namespace Feel
{
template <typename T, typename SizeT>
Backend<T,SizeT>::Backend( worldcomm_ptr_t const& worldComm )
    :
    super(worldComm),
#if defined( FEELPP_HAS_PETSC_H )
    M_backend    ( BACKEND_PETSC ),
#endif
    M_prefix( "" ),
    M_nlsolver(),
    M_rtolerance( 1e-8 ),
    M_dtolerance( 1e5 ),
    M_atolerance( 1e-50 ),
    M_rtoleranceSNES( 1e-8 ),
    M_stoleranceSNES( 1e-8 ),
    M_atoleranceSNES( 1e-50 ),
    M_rtoleranceKSPinSNES( 1e-5 ),
    M_reuse_prec( false ),
    M_reuse_jac( false ),
    M_reusePrecIsBuild(false),
    M_reusePrecRebuildAtFirstNewtonStep(false),
    M_reuseJacIsBuild(false),
    M_reuseJacRebuildAtFirstNewtonStep(true),
    M_transpose( false ),
    M_maxitKSP( 1000 ),
    M_maxitKSPinSNES( M_maxitKSP ),
    M_maxitSNES( 50 ),
    M_maxitKSPReuse( M_maxitKSP ),
    M_maxitKSPinSNESReuse( M_maxitKSPinSNES ),
    M_maxitSNESReuse( M_maxitSNES ),
    M_export( "" ),
    M_ksp( "gmres" ),
    M_pc( "lu" ),
    M_fieldSplit( "additive" ),
#if FEELPP_HAS_PETSC
    M_pcFactorMatSolverPackage( "petsc" ),
#endif
    M_constant_null_space( false ),
    M_showKSPMonitor( false ),
    M_showKSPConvergedReason( false ),
    M_verbose( false )
{
    if ( this->worldComm().globalSize() > 1 )
        M_pc = "gasm";
}

template <typename T, typename SizeT>
Backend<T,SizeT>::Backend( po::variables_map const& vm, std::string const& prefix, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    M_vm( vm ),
    M_prefix( prefix ),
    M_nlsolver( solvernonlinear_type::build( prefix, worldComm, vm ) ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_rtolerance( doption(_name="ksp-rtol",_prefix=prefix,_vm=vm) ),
    M_dtolerance( doption(_name="ksp-dtol",_prefix=prefix,_vm=vm) ),
    M_atolerance( doption(_name="ksp-atol",_prefix=prefix,_vm=vm) ),
    M_rtoleranceSNES( doption(_name="snes-rtol",_prefix=prefix,_vm=vm) ),
    M_stoleranceSNES( doption(_name="snes-stol",_prefix=prefix,_vm=vm) ),
    M_atoleranceSNES( doption(_name="snes-atol",_prefix=prefix,_vm=vm) ),
    M_rtoleranceKSPinSNES( doption(_name="snes-ksp-rtol",_prefix=prefix,_vm=vm) ),
    M_reuse_prec( boption(_name="reuse-prec",_prefix=prefix,_vm=vm) ),
    M_reuse_jac( boption(_name="reuse-jac",_prefix=prefix,_vm=vm) ),
    M_reusePrecIsBuild( false ) ,
    M_reusePrecRebuildAtFirstNewtonStep( boption(_name="reuse-prec.rebuild-at-first-newton-step",_prefix=prefix,_vm=vm) ),
    M_reuseJacIsBuild( false ) ,
    M_reuseJacRebuildAtFirstNewtonStep( boption(_name="reuse-jac.rebuild-at-first-newton-step",_prefix=prefix,_vm=vm) ),
    M_transpose( false ),
    M_maxitKSP( ioption(_name="ksp-maxit",_prefix=prefix,_vm=vm) ),
    M_maxitKSPinSNES( ioption(_name="snes-ksp-maxit",_prefix=prefix,_vm=vm) ),
    M_maxitSNES( ioption(_name="snes-maxit",_prefix=prefix,_vm=vm) ),
    M_maxitKSPReuse( (vm.count(prefixvm( prefix,"ksp-maxit-reuse")))?  ioption(_name="ksp-maxit-reuse",_prefix=prefix,_vm=vm) : M_maxitKSP ),
    M_maxitKSPinSNESReuse( (vm.count(prefixvm( prefix,"snes-ksp-maxit-reuse")))? ioption(_name="snes-ksp-maxit-reuse",_prefix=prefix,_vm=vm) : M_maxitKSPinSNES ),
    M_maxitSNESReuse( (vm.count(prefixvm( prefix,"snes-maxit-reuse")))? ioption(_name="snes-maxit-reuse",_prefix=prefix,_vm=vm) : M_maxitSNES ),
    M_export( soption(_name="export-matlab",_prefix=prefix,_vm=vm) ),
    M_ksp( soption(_name="ksp-type",_prefix=prefix,_vm=vm) ),
    M_pc( soption(_name="pc-type",_prefix=prefix,_vm=vm) ),
    M_fieldSplit( soption(_name="fieldsplit-type",_prefix=prefix,_vm=vm) ),
    M_pcFactorMatSolverPackage( soption(_name="pc-factor-mat-solver-package-type",_prefix=prefix,_vm=vm) ),
    M_constant_null_space( boption(_name="constant-null-space",_prefix=prefix,_vm=vm) ),
    M_showKSPMonitor( boption(_name="ksp-monitor",_prefix=prefix,_vm=vm) ),
    M_showKSPConvergedReason( vm.count(prefixvm( prefix,"ksp-converged-reason" )) ),
    M_verbose( boption(_name="backend.verbose",_prefix=prefix,_vm=vm) )
{
}
template <typename T, typename SizeT>
Backend<T,SizeT>::~Backend()
{
    this->clear();
}
template <typename T, typename SizeT>
void
Backend<T,SizeT>::clear()
{
    if ( M_preconditioner )
        M_preconditioner->clear();
    LOG(INFO) << "Sending delete signal to all observers... " << M_prefix << "\n";
    this->sendDeleteSignal();
    //this->clear ();
}
template <typename T, typename SizeT>
typename Backend<T,SizeT>::backend_ptrtype
Backend<T,SizeT>::build( BackendType bt, worldcomm_ptr_t const&worldComm )
{
    // Build the appropriate solver
    switch ( bt )
    {
    case BACKEND_EIGEN:
    {
        return backend_ptrtype( new BackendEigen<value_type>( worldComm ) );
    }
    break;
    case BACKEND_EIGEN_DENSE:
    {
        return backend_ptrtype( new BackendEigen<value_type,1>( worldComm ) );
    }
    break;

#if defined ( FEELPP_HAS_PETSC_H )

    case BACKEND_PETSC:
    {
        auto b = backend_ptrtype( new BackendPetsc<value_type>( worldComm ) );
        b->attachPreconditioner();
        return b;
    }
    break;
#endif

#if defined ( FEELPP_HAS_TRILINOS_EPETRA )

    case BACKEND_TRILINOS:
    {
        return backend_ptrtype( new BackendTrilinos( worldComm ) );
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

template <>
typename Backend<std::complex<double>>::backend_ptrtype
Backend<std::complex<double>>::build( BackendType bt, worldcomm_ptr_t const& worldComm )
{
    // Build the appropriate solver
    switch ( bt )
    {
    case BACKEND_EIGEN:
    {
        return backend_ptrtype( new BackendEigen<value_type>( worldComm ) );
    }
    break;
    case BACKEND_EIGEN_DENSE:
    {
        return backend_ptrtype( new BackendEigen<value_type,1>( worldComm ) );
    }
    break;

    default:
        std::cerr << "ERROR:  Unrecognized backend type package: "
                  << bt
                  << std::endl;
        throw std::invalid_argument( "invalid backend type" );
    }

    return backend_ptrtype();
}

template <typename T, typename SizeT>
typename Backend<T,SizeT>::backend_ptrtype
Backend<T,SizeT>::build( po::variables_map const& vm, std::string const& prefix, worldcomm_ptr_t const& worldComm )
{
    std::string kind = soption( _name="backend" );
    return build( kind, prefix, worldComm );
}
template <typename T, typename SizeT>
typename Backend<T,SizeT>::backend_ptrtype
Backend<T,SizeT>::build( std::string const& kind, std::string const& prefix, worldcomm_ptr_t const& worldComm, po::variables_map const& vm )
{
    if ( kind == "eigen")
        return backend_ptrtype( new BackendEigen<value_type,0,SizeT>( vm, prefix, worldComm ) );
    if ( kind == "eigen_dense")
        return backend_ptrtype( new BackendEigen<value_type,1,SizeT>( vm, prefix, worldComm ) );
#if defined ( FEELPP_HAS_PETSC_H )
    if ( kind == "petsc")
    {
        auto b = backend_ptrtype( new BackendPetsc<value_type,SizeT>( vm, prefix, worldComm ) );
        b->attachPreconditioner();
        return b;
    }
#else
    if ( kind == "petsc")
        LOG(FATAL) << "Backend 'petsc' not available";
#endif
    // should never happen
    return backend_ptrtype();
}

template <>
typename Backend<std::complex<double>>::backend_ptrtype
Backend<std::complex<double>>::build( std::string const& kind, std::string const& prefix, worldcomm_ptr_t const& worldComm, po::variables_map const& vm )
{
    if ( kind == "eigen")
        return backend_ptrtype( new BackendEigen<value_type,0,size_type>( vm, prefix, worldComm ) );
    if ( kind == "eigen_dense")
        return backend_ptrtype( new BackendEigen<value_type,1,size_type>( vm, prefix, worldComm ) );
    // should never happen
    return backend_ptrtype();
}


template <typename T, typename SizeT>
typename Backend<T,SizeT>::backend_ptrtype
Backend<T,SizeT>::build( BackendType bt, std::string const& prefix, worldcomm_ptr_t const& worldComm )
{
    return build( enumToKind( bt ), prefix, worldComm );
}

template<>
inline void
Backend<double>::attachPreconditioner()
{
    auto p = Feel::preconditioner( _prefix=this->prefix(),
                                   _pc=this->pcEnumType(),
                                   _backend=this->shared_from_this(),
                                   _pcfactormatsolverpackage=this->matSolverPackageEnumType() );
    if ( M_preconditioner && M_preconditioner != p )
        M_preconditioner->clear();
    M_preconditioner = p;
}

template <typename T, typename SizeT>
typename Backend<T,SizeT>::solve_return_type
Backend<T,SizeT>::solve( sparse_matrix_ptrtype const& A,
                   sparse_matrix_ptrtype const& P,
                   vector_ptrtype& x,
                   vector_ptrtype const& b,
                   bool reusePC )
{
    M_reusePC = reusePC;

    MatrixStructure matStructInitial = this->precMatrixStructure();
    size_type maxitKSPInitial = M_maxitKSP;
    vector_ptrtype x_save;

    if ( !M_reusePC )
    {
        //reset();
    }
    else
    {
        // save current solution in case of failure
        x->close();
        x_save = this->newVector(x->mapPtr());
        *x_save=*x;
        this->setPrecMatrixStructure( SAME_PRECONDITIONER );
        M_maxitKSP = M_maxitKSPReuse;
    }

    //start();


    //std::cout << "backend: " << this->precMatrixStructure() << "\n";
    boost::tie( M_converged, M_iteration, M_residual ) = this->solve( A, P, x, b );
    //stop();
    M_reuseFailed = reusePC && (!M_converged );

    if ( M_reuseFailed )
    {
        this->comm().globalComm().barrier();
        //reset();
        //start();

        // reset to initial solution
        x_save->close();
        *x=*x_save;

        this->setPrecMatrixStructure( matStructInitial );//DIFFERENT_NONZERO_PATTERN,SAME_NONZERO_PATTERN
        M_maxitKSP = maxitKSPInitial;

        if (this->comm().isMasterRank() )
            std::cout << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";
        LOG(INFO) << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";

        boost::tie( M_converged, M_iteration, M_residual ) = this->solve( A, P, x, b );

        //if ( !M_converged ) throw std::logic_error( "solver failed to converge" );
        if ( !M_converged ) std::cerr<< "Backend " << M_prefix << " : linear solver failed to converge" << std::endl;

        //stop();
    }

    this->setPrecMatrixStructure( matStructInitial );
    M_maxitKSP = maxitKSPInitial;

    return boost::make_tuple( M_converged, M_iteration, M_residual );
}
template <typename T, typename SizeT>
typename Backend<T,SizeT>::nl_solve_return_type
Backend<T,SizeT>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its,
                     bool reusePC, bool reuseJac )
{
    MatrixStructure matStructInitial = this->precMatrixStructure();

    M_nlsolver->setPreconditionerType( this->pcEnumType() );
    M_nlsolver->setKspSolverType( this->kspEnumType() );
    M_nlsolver->setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );

    M_nlsolver->setNbItMax( this->maxIterationsSNES() );
    M_nlsolver->setRelativeResidualTol( this->rToleranceSNES() );
    M_nlsolver->setAbsoluteResidualTol( this->aToleranceSNES() );
    M_nlsolver->setAbsoluteSolutionTol( this->sToleranceSNES() );
    M_nlsolver->setRtoleranceKSP( this->rtoleranceKSPinSNES() );
    M_nlsolver->setAtoleranceKSP( this->aTolerance() );
    M_nlsolver->setDtoleranceKSP( this->dTolerance() );
    M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNES() );

    M_nlsolver->attachNullSpace( this->M_nullSpace );
    M_nlsolver->attachNearNullSpace( this->M_nearNullSpace );

    M_nlsolver->init();

    //vector_ptrtype x_save = x->clone();
    vector_ptrtype x_save;

    //std::cout << "[nlSolve] reusepc:" << reusePC << std::endl;
    if ( reusePC || reuseJac )
    {
        // save current solution in case of failure
        x->close();
        x_save = this->newVector(x->mapPtr());
        *x_save=*x;

        // configure reusePC,reuseJac in non linear solver
        int typeReusePrec = 1,typeReuseJac = 1 ;
#if FEELPP_HAS_PETSC
#if PETSC_VERSION_LESS_THAN(3,5,0)
        // if first time or rebuild prec at first newton step, need to get matStructInitial
        if ( reusePC && (!M_reusePrecIsBuild || M_reusePrecRebuildAtFirstNewtonStep) )
            {
                M_nlsolver->setPrecMatrixStructure( matStructInitial );
                M_reusePrecIsBuild=true;
            }
        else if ( reusePC ) M_nlsolver->setPrecMatrixStructure( SAME_PRECONDITIONER );

        if ( reusePC ) typeReusePrec = -1;
#else
        if ( reusePC )
        {
            M_nlsolver->setPrecMatrixStructure( SAME_PRECONDITIONER );
            typeReusePrec = (M_reusePrecRebuildAtFirstNewtonStep)? -2 : -1;
        }
#endif
#endif
        if ( reuseJac ) typeReuseJac = (M_reuseJacRebuildAtFirstNewtonStep)? -2 : -1;

        M_nlsolver->setReuse( typeReuseJac, typeReusePrec );

        // special tolerance in reuse mode
        M_nlsolver->setNbItMax( this->maxIterationsSNESReuse() );
        M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNESReuse() );

        // compute cst jacobian when jacobian is never rebuilt after!
        if ( reuseJac &&  (!M_reuseJacIsBuild && !M_reuseJacRebuildAtFirstNewtonStep) ) { this->nlSolver()->jacobian( x, A );A->close();M_reuseJacIsBuild=true;}
    }
    else
    {
        M_nlsolver->setReuse( 1, 1 );
    }

    auto ret = M_nlsolver->solve( A, x, b, tol, its );

    //std::cout << "[nlSolve] ret.first " << ret.first <<std::endl;
    if ( !ret.isConverged() && ( reusePC || reuseJac ) )
    {
        if ( this->comm().isMasterRank() )
            std::cout << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";
        LOG(INFO) << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";

        // reset to initial solution
        x_save->close();
        *x=*x_save;
        x->close();

        //M_nlsolver->init();
        //M_nlsolver->setPreconditionerType( this->pcEnumType() );
        //M_nlsolver->setKspSolverType( this->kspEnumType() );
        M_nlsolver->setPrecMatrixStructure( matStructInitial/*SAME_NONZERO_PATTERN*/ );
        M_nlsolver->setNbItMax( this->maxIterationsSNES() );
        M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNES() );

        //M_nlsolver->setReuse( 1, -2 );
        M_nlsolver->setReuse( 1, 1 );

        this->nlSolver()->jacobian( x, A );
        A->close();
        // call solver which must execute with success
        auto ret2 = M_nlsolver->solve( A, x, b, tol, its );

        if ( !ret2.isConverged() )
        {
            LOG(ERROR) << "\n[backend] non-linear solver fail";
            LOG(ERROR) << "Backend " << M_prefix << " : non-linear solver failed to converge" << std::endl;
        }

        return ret2;
    }
    else if ( !ret.isConverged() )
    {
        if ( this->worldComm().isMasterRank() )
        {
            LOG(ERROR) << "\n[backend] non-linear solver fail";
            LOG(ERROR) << "Backend " << M_prefix << " : non-linear solver failed to converge" << std::endl;
        }
    }

    this->setPrecMatrixStructure( matStructInitial );

    return ret;
}
template <typename T, typename SizeT>
typename Backend<T,SizeT>::nl_solve_return_type
Backend<T,SizeT>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its )
{

    M_nlsolver->setPreconditionerType( this->pcEnumType() );
    M_nlsolver->setKspSolverType( this->kspEnumType() );
    M_nlsolver->setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );
    M_nlsolver->setNbItMax( this->maxIterationsSNES() );
    M_nlsolver->setRelativeResidualTol( this->rToleranceSNES() );
    M_nlsolver->setAbsoluteResidualTol( this->aToleranceSNES() );
    M_nlsolver->setAbsoluteSolutionTol( this->sToleranceSNES() );
    M_nlsolver->setRtoleranceKSP( this->rtoleranceKSPinSNES() );
    M_nlsolver->setAtoleranceKSP( this->aTolerance() );
    M_nlsolver->setDtoleranceKSP( this->dTolerance() );
    M_nlsolver->setMaxitKSP( this->maxIterationsKSP() );

    M_nlsolver->init();
    M_nlsolver->setReuse( 1, 1 );

    auto ret = M_nlsolver->solve( A, x, b, tol, its );

    if ( !ret.isConverged() )
    {
        LOG(ERROR) << "\n[backend] non-linear solver fail";
        LOG(ERROR) << "Backend " << M_prefix << " : non-linear solver failed to converge" << std::endl;
    }

    return ret;
}
template <typename T, typename SizeT>
int
Backend<T,SizeT>::PtAP( sparse_matrix_ptrtype const& A,
                  sparse_matrix_ptrtype const& P,
                  sparse_matrix_ptrtype & C ) const
{
    LOG(WARNING) << "PtAP not implemented in base class. You need to implement the procedure in the current backend.";
    return 0;
}
template <typename T, typename SizeT>
int
Backend<T,SizeT>::PAPt( sparse_matrix_ptrtype const& A,
                  sparse_matrix_ptrtype const& P,
                  sparse_matrix_ptrtype & C ) const
{
    LOG(WARNING) << "PAPt not implemented in base class. You need to implement the procedure in the current backend.";
    return 0;
}

template <typename T, typename SizeT>
typename Backend<T,SizeT>::value_type
Backend<T,SizeT>::dot( vector_type const& x, vector_type const& y ) const
{
    return x.dot( y );
}

template <>
typename Backend<std::complex<double>>::value_type
Backend<std::complex<double>>::dot( vector_type const& x, vector_type const& y ) const
{
    value_type localres = 0;

    for ( size_type i = 0; i < x.map().nLocalDofWithoutGhost(); ++i )
    {
        localres += std::conj(x( i ))*y( i );
    }

    value_type globalres=localres;
    mpi::all_reduce( this->worldComm().globalComm(), localres, globalres, std::plus<value_type>() );
    return globalres;
}

template <typename T, typename SizeT>
void
Backend<T,SizeT>::start()
{
    M_timer.start();
}

template <typename T, typename SizeT>
void
Backend<T,SizeT>::stop()
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
            M_maxitKSP = std::max(size_type(10),std::min( M_maxitKSP, ( size_type )( 1.5*solveIter + 10.5 ) ));
    }

    else
    {
        double nextSolveIter;

        if ( M_nUsePC == 2 )
        {
            M_totalSolveIter = solveIter*( 1.0+M_firstSolveTime/solveTime );
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
            M_maxitKSP = std::max( size_type(10),std::min( M_maxitKSP, ( size_type )( M_totalSolveIter/M_nUsePC + 0.5 ) ));
        }
    }
}

template <typename T, typename SizeT>
void
Backend<T,SizeT>::reset()
{
    M_reusePC = false;
    M_totalSolveIter = 0.0;
    M_nUsePC = 0;
    //M_backend->set_maxiter( M_maxitKSP );

}

template <typename T, typename SizeT>
SolverType
Backend<T,SizeT>::kspEnumType() const
{
    return kspTypeConvertStrToEnum( this->kspType() );
}

template <typename T, typename SizeT>
SolverNonLinearType
Backend<T,SizeT>::snesEnumType() const
{
    return snesTypeConvertStrToEnum( this->snesType() );
}

template <typename T, typename SizeT>
PreconditionerType
Backend<T,SizeT>::pcEnumType() const
{
    return pcTypeConvertStrToEnum( this->pcType() );
}

template <typename T, typename SizeT>
FieldSplitType
Backend<T,SizeT>::fieldSplitEnumType() const
{
    return fieldsplitTypeConvertStrToEnum( this->fieldsplitType() );
}

template <typename T, typename SizeT>
MatSolverPackageType
Backend<T,SizeT>::matSolverPackageEnumType() const
{
    return matSolverPackageConvertStrToEnum( this->pcFactorMatSolverPackageType() );
}


#if 0
void updateBackendPreconditionerOptions( po::options_description & _options, std::string const& prefix, std::string const& sub = "",
                                         std::string pcType = "lu", bool useDefaultValue=true )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    _options.add_options()
        ( prefixvm( prefix,pcctx+"pc-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( pcType ):Feel::po::value<std::string>(),
          "type of preconditioners (lu, ilut, ilutp, diag, id,...)" )
        ( prefixvm( prefix,pcctx+"pc-view" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "display preconditioner information" )
        ( prefixvm( prefix,pcctx+"pc-use-config-default-petsc" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "configure pc with default petsc options" )
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
        ( prefixvm( prefix,pcctx+"pc-factor-mat-solver-package-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "mumps" ):Feel::po::value<std::string>(),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#else
        ( prefixvm( prefix,pcctx+"pc-factor-mat-solver-package-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "petsc" ):Feel::po::value<std::string>(),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif
        ( prefixvm( prefix,pcctx+"ilu-threshold" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1e-3 ):Feel::po::value<double>(),
          "threshold value for preconditioners" )
        ( prefixvm( prefix,pcctx+"ilu-fillin" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 2 ):Feel::po::value<int>(),
          "fill-in level value for preconditioners" )
        ( prefixvm( prefix,pcctx+"pc-factor-levels" ).c_str(),
        (useDefaultValue)?Feel::po::value<int>()->default_value( 3 ):Feel::po::value<int>(),
        "Sets the number of levels of fill to use for ilu" )
        ( prefixvm( prefix,pcctx+"pc-factor-fill" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 6 ):Feel::po::value<double>(),
          "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )

        ( prefixvm( prefix,pcctx+"pc-sor-omega" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1. ):Feel::po::value<double>(),
          "Sets the SOR relaxation coefficient, omega" )
        ( prefixvm( prefix,pcctx+"pc-sor-lits" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of local iterations, smoothings over just variables on processor" )
        ( prefixvm( prefix,pcctx+"pc-sor-its" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of parallel iterations to use; each parallel iteration has lits local iterations" )
        ( prefixvm( prefix,pcctx+"pc-sor-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "local_symmetric" ):Feel::po::value<std::string>(),
          "(symmetric,forward,backward,local_symmetric,local_forward,local_backward) Sets the SOR preconditioner to use symmetric (SSOR), backward, or forward relaxation. The local variants perform SOR on each processor" )
        ;

#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
        for ( int icntl=1 ; icntl<= 33 ; ++icntl )
        {
            std::string mumpsOption = (boost::format("pc-factor-mumps.icntl-%1%")%icntl ).str();
            _options.add_options()
                ( prefixvm( prefix,pcctx+mumpsOption ).c_str(),
                  Feel::po::value<int>(),"configure mumps factorisation (see mumps ICNTL documentation)" );
        }
#endif


    // add gasm and asm (can not used as sub preconditioner)
    if ( sub.empty() )
        _options.add_options()
            ( prefixvm( prefix,"pc-gasm-type" ).c_str(), (useDefaultValue)?Feel::po::value<std::string>()->default_value( "restrict" ):Feel::po::value<std::string>(),
              "type of gasm (basic, restrict, interpolate, none)" )
            ( prefixvm( prefix,"pc-gasm-overlap" ).c_str(), (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
              "number of overlap levels" )
            ( prefixvm( prefix,"pc-asm-type" ).c_str(), (useDefaultValue)?Feel::po::value<std::string>()->default_value( "restrict" ):Feel::po::value<std::string>(),
              "type of asm (basic, restrict, interpolate, none)" )
            ( prefixvm( prefix,"pc-asm-overlap" ).c_str(), (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
              "number of overlap levels" )
            ;
}
#endif
void updateBackendKSPOptions( po::options_description & _options, std::string const& prefix, std::string const& sub = "",
                              std::string const& kspType = "gmres",double rtol = 1e-8, size_type maxit=1000, bool useDefaultValue=true  )
{
    std::string kspctx = (sub.empty())? "" : sub+"-";
    _options.add_options()
        ( prefixvm( prefix,kspctx+"ksp-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( kspType ):Feel::po::value<std::string>(),
          "cg, bicgstab, gmres,preonly,..." )
        ( prefixvm( prefix,kspctx+"ksp-view" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "Prints the KSP data structure" )
        ( prefixvm( prefix,kspctx+"ksp-monitor" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "monitor ksp" )
        ( prefixvm( prefix,kspctx+"ksp-converged-reason" ).c_str() , "converged reason ksp" )
        ( prefixvm( prefix,kspctx+"ksp-verbose" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 0 ):Feel::po::value<int>(),
          "(=0,1,2) print solver iterations" )
        ( prefixvm( prefix,kspctx+"ksp-rtol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( rtol ):Feel::po::value<double>(),
          "relative tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-atol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1e-50 ):Feel::po::value<double>(),
          "absolute tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-dtol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1e5 ):Feel::po::value<double>(),
          "divergence tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-maxit" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( maxit ):Feel::po::value<int>(),
          "maximum number of iterations" )
        ( prefixvm( prefix,kspctx+"ksp-maxit-reuse" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>():Feel::po::value<int>(),
          "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,kspctx+"constant-null-space" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( 0 ):Feel::po::value<bool>(),
          "set the null space to be the constant values" )
#if FEELPP_HAS_PETSC
        ( prefixvm( prefix,kspctx+"ksp-use-config-default-petsc" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "configure ksp with default petsc options" )
#endif

        ( prefixvm( prefix,kspctx+"ksp-use-initial-guess-nonzero" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "tells the iterative solver that the initial guess is nonzero" )

        // Default value : "" let the default behavior
#if FEELPP_HAS_PETSC
        ( prefixvm( prefix,kspctx+"ksp-norm-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "default" ),
          "Sets the norm that is used for convergence testing. Leaving this empty will use default petsc behavior. default, none, preconditioned, unpreconditioned, natural" )
#endif

        ( prefixvm( prefix,kspctx+"gmres-restart" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 30 ):Feel::po::value<int>(),
          "number of iterations before solver restarts (gmres)" )
        ( prefixvm( prefix,kspctx+"fgmres-restart" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 30 ):Feel::po::value<int>(),
          "number of iterations before solver restarts (fgmres)" )
        ( prefixvm( prefix,kspctx+"gcr-restart" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 30 ):Feel::po::value<int>(),
          "number of iterations before solver restarts (gcr)" )
        ;
}
#if 0
void updateBackendMGPreconditionerOptions( po::options_description & _options, std::string const& prefix, std::string const& sub = "" )
{
    std::string pcctx = (sub.empty())? "pc-" : sub+"-pc-";
    // multigrid options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"mg-levels" ).c_str(), Feel::po::value<int>()->default_value( 10/*2*/ ), "number of levels including finest" )
        ( prefixvm( prefix,pcctx+"mg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "multiplicative" ), "Determines the form of multigrid to use: multiplicative, additive, full, kaskade " )
        ( prefixvm( prefix,pcctx+"mg-smoothdown" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of smoothing steps before applying restriction operator" )
        ;
    // ml options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"ml-reuse-interpolation" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Reuse the interpolation operators when possible (cheaper, weaker when matrix entries change a lot)" )
        ( prefixvm( prefix,pcctx+"ml-keep-agg-info" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Allows the preconditioner to be reused, or auxiliary matrices to be generated" )
        ( prefixvm( prefix,pcctx+"ml-reusable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Store intermedaiate data structures so that the multilevel hierarchy is reusable" )
        ( prefixvm( prefix,pcctx+"ml-old-hierarchy" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Use old routine to generate hierarchy" )
        ;
    // gamg options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"gamg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "agg" ), "type of generalized algebraic multigrid : agg, geo " )
        ( prefixvm( prefix,pcctx+"gamg-proc-eq-lim" ).c_str(), Feel::po::value<int>()->default_value( 50 ), "number of equations to aim for on coarse grids via processor reduction" )
        ( prefixvm( prefix,pcctx+"gamg-coarse-eq-lim" ).c_str(), Feel::po::value<int>()->default_value( 800 ), "max number of equations on coarse grids" )
        ( prefixvm( prefix,pcctx+"gamg-threshold" ).c_str(), Feel::po::value<double>()->default_value( 0. ), "relative threshold to use for dropping edges in aggregation graph" )
        ;
    // coarse ksp/pc
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";
    std::string prefixMGCoarse = ( boost::format( "%1%%2%coarse" ) %prefixvm(prefix,"") %mgctx ).str();
    updateBackendPreconditionerOptions( _options, prefixMGCoarse, "", "redundant"/*:"lu"*/ );
    updateBackendPreconditionerOptions( _options, prefixMGCoarse, "sub", "lu" );
    updateBackendKSPOptions( _options,prefixMGCoarse,"", "preonly",1e-5,50);
    updateBackendKSPOptions( _options,prefixMGCoarse,"sub", "preonly",1e-5,50);
    // all levels ksp/pc (not including coarse level) with default values
    std::string prefixMGLevelsGeneric = prefixvm( prefix, mgctx+"levels" );
    updateBackendPreconditionerOptions( _options, prefixMGLevelsGeneric,    "", "sor" );
    updateBackendPreconditionerOptions( _options, prefixMGLevelsGeneric, "sub", "lu"  ); // gasm/asm
    updateBackendKSPOptions( _options,prefixMGLevelsGeneric,   "", "richardson", 1e-5, 2  );
    updateBackendKSPOptions( _options,prefixMGLevelsGeneric,"sub",    "preonly", 1e-5, 50 ); // gasm/asm
    // fine level
    std::string prefixMGFineLevel = prefixvm( prefix, mgctx+"fine-level" );
    updateBackendPreconditionerOptions( _options, prefixMGFineLevel,    "", "sor", false );
    updateBackendPreconditionerOptions( _options, prefixMGFineLevel, "sub", "lu", false  ); // gasm/asm
    updateBackendKSPOptions( _options, prefixMGFineLevel,    "", "richardson", 1e-5, 2, false  );
    updateBackendKSPOptions( _options, prefixMGFineLevel, "sub",    "preonly", 1e-5, 50, false ); // gasm/asm
    // each levels 1 to 5 can be control separately
    for ( uint16_type i=1; i<6; ++i )
    {
        std::string prefixMGLevels = ( boost::format( "%1%%2%levels%3%" ) %prefixvm(prefix,"") %mgctx %i ).str();
        updateBackendPreconditionerOptions( _options, prefixMGLevels,    "", "sor", false );
        updateBackendPreconditionerOptions( _options, prefixMGLevels, "sub",  "lu", false ); // gasm/asm
        updateBackendKSPOptions( _options,prefixMGLevels,   "", "richardson", 1e-5,  2, false  );
        updateBackendKSPOptions( _options,prefixMGLevels,"sub",    "preonly", 1e-5, 50, false ); // gasm/asm
    }
}

void updateBackendFieldSplitPreconditionerOptions( po::options_description & _options, std::string const& prefix, std::string const& sub = "" )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    // field split options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"fieldsplit-type" ).c_str(), Feel::po::value<std::string>()->default_value( "additive" ), "type of fieldsplit (additive, multiplicative, symmetric-multiplicative, schur)" )
        ( prefixvm( prefix,pcctx+"fieldsplit-fields" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "fields definition (ex: --fieldsplit-fields=0->(0,2),1->(1)" )
        ;
    // schur complement options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"fieldsplit-schur-fact-type" ).c_str(), Feel::po::value<std::string>()->default_value( "full" ), "type of schur factorization (diag, lower, upper, full)" )
        ( prefixvm( prefix,pcctx+"fieldsplit-schur-precondition" ).c_str(), Feel::po::value<std::string>()->default_value( "a11" ), "self,user,a11" )
        ;

    // inner solver (A^{-1}) of schur complement : S = C-B A^{-1} B^T
    std::string prefixSchurInnerSolver = prefixvm( prefix,pcctx+"fieldsplit-schur-inner-solver" );
    _options.add_options()
        ( prefixvm( prefixSchurInnerSolver,"use-outer-solver" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use-outer-solver" )
        ;
    updateBackendPreconditionerOptions( _options, prefixSchurInnerSolver ,"", "jacobi" );
    updateBackendPreconditionerOptions( _options, prefixSchurInnerSolver, "sub", "lu" ); // gasm
    updateBackendKSPOptions( _options, prefixSchurInnerSolver,   "", "preonly", 1e-5,  10, true  ); // preonly or gmres??
    updateBackendKSPOptions( _options, prefixSchurInnerSolver, "sub", "preonly", 1e-5, 50 ); // ksp with gasm

    // solver (A^{-1}) used in upper schur preconditioning
    std::string prefixSchurUpperSolver = prefixvm( prefix,pcctx+"fieldsplit-schur-upper-solver" );
    _options.add_options()
        ( prefixvm( prefixSchurUpperSolver,"use-outer-solver" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use-outer-solver" )
        ;
    updateBackendPreconditionerOptions( _options, prefixSchurUpperSolver ,"", "jacobi" );
    updateBackendKSPOptions( _options, prefixSchurUpperSolver,   "", "preonly", 1e-5,  10, true  ); // preonly or gmres??

    // ksp/pc options for each split
    for ( uint16_type i=0; i<5; ++i )
    {
        std::string prefixfieldsplit = ( boost::format( "%1%%2%fieldsplit-%3%" ) %prefixvm( prefix,"" ) %pcctx %i ).str();
        updateBackendKSPOptions( _options, prefixfieldsplit, "", "preonly", 1e-5, 50 );
        updateBackendKSPOptions( _options, prefixfieldsplit, "sub", "preonly", 1e-5, 50 ); // ksp with gasm
        updateBackendPreconditionerOptions( _options, prefixfieldsplit, "", (i==0)?"lu":"none" ); // lu/ilu
        updateBackendPreconditionerOptions( _options, prefixfieldsplit, "sub", "lu" ); // gasm
        updateBackendMGPreconditionerOptions( _options, prefixfieldsplit ,"" ); // multigrid
    }

    // specific case : lsc preconditioner ( only with split 1 )
    std::string prefixfieldsplitLSC = prefixvm( prefixvm( prefix,pcctx+"fieldsplit-1" ), "lsc" );
    _options.add_options()
        ( prefixvm( prefixfieldsplitLSC,"scale-diag" ).c_str(), Feel::po::value<bool>()->default_value( false ), "scale diag" )
        ;
    updateBackendPreconditionerOptions( _options, prefixfieldsplitLSC ,"" );
    updateBackendPreconditionerOptions( _options, prefixfieldsplitLSC ,"sub" );
    updateBackendMGPreconditionerOptions( _options, prefixfieldsplitLSC ,"" ); // multigrid
    updateBackendKSPOptions( _options, prefixfieldsplitLSC, "", "preonly", 1e-5, 50 ); // lsc+(lu/ilu)
    updateBackendKSPOptions( _options, prefixfieldsplitLSC, "sub", "preonly", 1e-5, 50 ); // lsc+gasm

}
#endif
/**
 * \return the command lines options of the petsc backend
 */
po::options_description backend_options( std::string const& prefix )
{
    po::options_description _options( "Linear and NonLinear Solvers Backend " + prefix + " options" );
    _options.add_options()
        // solver options
        ( prefixvm( prefix,"backend" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ), "backend type: petsc, eigen, eigen_dense" )
        ( prefixvm( prefix,"backend.rebuild" ).c_str(), Feel::po::value<bool>()->default_value( false ), "rebuild the backend each time it is called" )
        ( prefixvm( prefix,"backend.rebuild_op" ).c_str(), Feel::po::value<bool>()->default_value( true ), "rebuild the backend associated to operators" )
        ( prefixvm( prefix,"backend.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "set the backend to be verbose" )

        ( prefixvm( prefix,"reuse-jac" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reuse jacobian" )
        ( prefixvm( prefix,"reuse-jac.rebuild-at-first-newton-step" ).c_str(), Feel::po::value<bool>()->default_value( true ), "rebuild jacobian at each Newton when reuse jacobian" )
        ( prefixvm( prefix,"reuse-prec" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reuse preconditioner" )
        ( prefixvm( prefix,"reuse-prec.rebuild-at-first-newton-step" ).c_str(), Feel::po::value<bool>()->default_value( false ), "rebuild preconditioner at each Newton when reuseprec" )

        ( prefixvm( prefix,"export-matlab" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "export matrix/vector to matlab, default empty string means no export, other string is used as prefix" )

        ( prefixvm( prefix,"snes-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Prints the SNES data structure" )
        ( prefixvm( prefix,"snes-type" ).c_str(), Feel::po::value<std::string>()->default_value( "ls" ), "Set the SNES solver" )
        ( prefixvm( prefix,"snes-line-search-type" ).c_str(), Feel::po::value<std::string>()->default_value( "bt" ), "Set the SNES line search solver" )
        ( prefixvm( prefix,"snes-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-8 ), "relative tolerance" )
        ( prefixvm( prefix,"snes-atol" ).c_str(), Feel::po::value<double>()->default_value( 1e-50 ), "absolute tolerance" )
        ( prefixvm( prefix,"snes-stol" ).c_str(), Feel::po::value<double>()->default_value( 1e-8 ), "step length tolerance" )
        ( prefixvm( prefix,"snes-maxit" ).c_str(), Feel::po::value<int>()->default_value( 50 ), "maximum number of iterations" )
        ( prefixvm( prefix,"snes-maxit-reuse" ).c_str(), Feel::po::value<int>(), "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,"snes-ksp-maxit" ).c_str(), Feel::po::value<int>()->default_value( 1000 ), "maximum number of iterations" )
        ( prefixvm( prefix,"snes-ksp-maxit-reuse" ).c_str(), Feel::po::value<int>(), "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,"snes-ksp-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-5 ), "relative tolerance" )
        ( prefixvm( prefix,"snes-monitor" ).c_str(), Feel::po::value<bool>()->default_value( false ) , "monitor snes" )
        ( prefixvm( prefix,"snes-converged-reason" ).c_str() , "converged reason snes" )
        ;

    updateBackendKSPOptions( _options, prefix, "", "gmres" );// ksp options
#if 0
    updateBackendPreconditionerOptions( _options, prefix, "", "lu" ); // pc options
    updateBackendKSPOptions( _options, prefix, "sub", "preonly" );// gasm/bjacobi + ksp
    updateBackendPreconditionerOptions( _options, prefix, "sub", "lu" ); // gasm/asm

    updateBackendMGPreconditionerOptions( _options, prefix ); // multigrid
    updateBackendMGPreconditionerOptions( _options, prefix, "sub" ); // (gasm/bjacobi)+multigrid

    updateBackendFieldSplitPreconditionerOptions( _options, prefix ); // fieldsplit
    updateBackendFieldSplitPreconditionerOptions( _options, prefix, "sub" ); // (gasm/bjacobi) + fieldsplit
    updateBackendFieldSplitPreconditionerOptions( _options, prefixvm( prefix,"fieldsplit-0" ) ); // fieldsplit + (fieldsplit in subsplit0)
#else
    _options.add_options()
        ( prefixvm( prefix,"pc-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "lu" ),
          "type of preconditioners (lu, ilut, ilutp, diag, id,...)" )
        ( prefixvm( prefix,"pc-view" ).c_str(),
          Feel::po::value<bool>()->default_value( false ),
          "display preconditioner information" )
#if FEELPP_HAS_PETSC
        ( prefixvm( prefix,"pc-use-config-default-petsc" ).c_str(),
          Feel::po::value<bool>()->default_value( false ),
          "configure pc with default petsc options" )
#endif
        ( prefixvm( prefix,"pc-factor-shift-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "none" ),
          "adds a particular type of quantity to the diagonal of the matrix during numerical factorization, thus the matrix has nonzero pivots (none, nonzero, positive_definite, inblocks)" )
#if defined(FEELPP_HAS_PETSC)
#if defined(PETSC_HAVE_MUMPS)
        ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "mumps" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#else
        ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "petsc" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif
#else
        ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(),
          Feel::po::value<std::string>()->default_value( "feelpp" /*"mumps"*/ ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif
        ( prefixvm( prefix,"fieldsplit-type" ).c_str(), Feel::po::value<std::string>()->default_value( "additive" ),
          "type of fieldsplit (additive, multiplicative, symmetric-multiplicative, schur)" )
        ( prefixvm( prefix,"fieldsplit-fields" ).c_str(), Feel::po::value<std::string>()->default_value( "" ),
          "fields definition (ex: --fieldsplit-fields=0->(0,2),1->(1)" )
        ( prefixvm( prefix,"fieldsplit-use-components" ).c_str(), Feel::po::value<bool>()->default_value( false ),"split also with components" )
        ;
#endif


    return _options;
}


/*
 * Explicit instantiations
 */
template class Backend<double,uint32_type>;
template class Backend<std::complex<double>,uint32_type>;


} // namespace Feel

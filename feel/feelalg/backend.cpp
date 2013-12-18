/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-12-23

  Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)

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
#include <feel/feelalg/backend.hpp>
// PETSc defines MatType which is used as a typename by eigen3 and it conflicts
// undef MatType here to ensure compilation since it is not needed to compile this file
#undef MatType
#include <feel/feelalg/backendeigen.hpp>
#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feelalg/backendtrilinos.hpp>


namespace Feel
{
template <typename T>
Backend<T>::Backend( WorldComm const& worldComm )
    :
    M_worldComm(worldComm),
#if defined( FEELPP_HAS_PETSC_H )
    M_backend    ( BACKEND_PETSC ),
#endif
    M_prefix( "" ),
    M_nlsolver(),
    M_rtolerance( 1e-13 ),
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
    M_snesType( "ls" ),
    M_pc( "lu" ),
    M_fieldSplit( "additive" ),
    M_pcFactorMatSolverPackage( "petsc" ),
    M_constant_null_space( false ),
    M_showKSPMonitor( false ), M_showSNESMonitor( false ),
    M_showKSPConvergedReason( false ), M_showSNESConvergedReason( false )
{
    if ( M_worldComm.globalSize() > 1 )
        M_pc = "gasm";
}

template <typename T>
Backend<T>::Backend( Backend const& backend )
    :
    M_worldComm( backend.M_worldComm ),
    M_backend( backend.M_backend ),
    M_prefix( backend.M_prefix ),
    M_nlsolver( backend.M_nlsolver ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_rtolerance( backend.M_rtolerance ),
    M_dtolerance( backend.M_dtolerance ),
    M_atolerance( backend.M_atolerance ),
    M_rtoleranceSNES( backend.M_rtoleranceSNES ),
    M_stoleranceSNES( backend.M_stoleranceSNES ),
    M_atoleranceSNES( backend.M_atoleranceSNES ),
    M_rtoleranceKSPinSNES( backend.M_rtoleranceKSPinSNES ),
    M_reuse_prec( backend.M_reuse_prec ),
    M_reuse_jac( backend.M_reuse_jac ),
    M_reusePrecIsBuild( backend.M_reusePrecIsBuild) ,
    M_reusePrecRebuildAtFirstNewtonStep( backend.M_reusePrecRebuildAtFirstNewtonStep ),
    M_reuseJacIsBuild( backend.M_reuseJacIsBuild) ,
    M_reuseJacRebuildAtFirstNewtonStep( backend.M_reuseJacRebuildAtFirstNewtonStep ),
    M_transpose( backend.M_transpose ),
    M_maxitKSP( backend.M_maxitKSP ),
    M_maxitKSPinSNES( backend.M_maxitKSPinSNES ),
    M_maxitSNES( backend.M_maxitSNES ),
    M_maxitKSPReuse( backend.M_maxitKSPReuse ),
    M_maxitKSPinSNESReuse( backend.M_maxitKSPinSNESReuse ),
    M_maxitSNESReuse( backend.M_maxitSNESReuse ),
    M_export( backend.M_export ),
    M_ksp( backend.M_ksp ),
    M_snesType( backend.M_snesType ),
    M_pc( backend.M_pc ),
    M_fieldSplit( backend.M_fieldSplit ),
    M_pcFactorMatSolverPackage( backend.M_pcFactorMatSolverPackage ),
    M_constant_null_space( backend.M_constant_null_space ),
    M_showKSPMonitor( backend.M_showKSPMonitor ),
    M_showSNESMonitor( backend.M_showSNESMonitor ),
    M_showKSPConvergedReason( backend.M_showKSPConvergedReason ),
    M_showSNESConvergedReason( backend.M_showSNESConvergedReason )
{
}
template <typename T>
Backend<T>::Backend( po::variables_map const& vm, std::string const& prefix, WorldComm const& worldComm )
    :
    M_worldComm( worldComm ),
    M_vm( vm ),
    M_prefix( prefix ),
    M_nlsolver( solvernonlinear_type::build( vm, prefix, worldComm ) ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_rtolerance( vm[prefixvm( prefix,"ksp-rtol" )].template as<double>() ),
    M_dtolerance( vm[prefixvm( prefix,"ksp-dtol" )].template as<double>() ),
    M_atolerance( vm[prefixvm( prefix,"ksp-atol" )].template as<double>() ),
    M_rtoleranceSNES( vm[prefixvm( prefix,"snes-rtol" )].template as<double>() ),
    M_stoleranceSNES( vm[prefixvm( prefix,"snes-stol" )].template as<double>() ),
    M_atoleranceSNES( vm[prefixvm( prefix,"snes-atol" )].template as<double>() ),
    M_rtoleranceKSPinSNES( vm[prefixvm( prefix,"snes-ksp-rtol" )].template as<double>() ),
    M_reuse_prec( vm[prefixvm( prefix,"reuse-prec" )].template as<bool>() ),
    M_reuse_jac(  vm[prefixvm( prefix,"reuse-jac" )].template as<bool>() ),
    M_reusePrecIsBuild( false ) ,
    M_reusePrecRebuildAtFirstNewtonStep( vm[prefixvm( prefix,"reuse-prec.rebuild-at-first-newton-step" )].template as<bool>() ),
    M_reuseJacIsBuild( false ) ,
    M_reuseJacRebuildAtFirstNewtonStep( vm[prefixvm( prefix,"reuse-jac.rebuild-at-first-newton-step" )].template as<bool>() ),
    M_transpose( false ),
    M_maxitKSP( vm[prefixvm( prefix,"ksp-maxit" )].template as<size_type>() ),
    M_maxitKSPinSNES( vm[prefixvm( prefix,"snes-ksp-maxit" )].template as<size_type>() ),
    M_maxitSNES( vm[prefixvm( prefix,"snes-maxit" )].template as<size_type>() ),
    M_maxitKSPReuse( (vm.count(prefixvm( prefix,"ksp-maxit-reuse")))? vm[prefixvm( prefix,"ksp-maxit-reuse" )].template as<size_type>() : M_maxitKSP ),
    M_maxitKSPinSNESReuse( (vm.count(prefixvm( prefix,"snes-ksp-maxit-reuse")))? vm[prefixvm( prefix,"snes-ksp-maxit-reuse" )].template as<size_type>() : M_maxitKSPinSNES ),
    M_maxitSNESReuse( (vm.count(prefixvm( prefix,"snes-maxit-reuse")))? vm[prefixvm( prefix,"snes-maxit-reuse" )].template as<size_type>() : M_maxitSNES ),
    M_export( vm[prefixvm( prefix,"export-matlab" )].template as<std::string>() ),
    M_ksp( vm[prefixvm( prefix,"ksp-type" )].template as<std::string>() ),
    M_snesType( vm[prefixvm( prefix,"snes-type" )].template as<std::string>() ),
    M_pc( vm[prefixvm( prefix,"pc-type" )].template as<std::string>() ),
    M_fieldSplit( vm[prefixvm( prefix,"fieldsplit-type" )].template as<std::string>() ),
    M_pcFactorMatSolverPackage( vm[prefixvm( prefix,"pc-factor-mat-solver-package-type" )].template as<std::string>() ),
    M_constant_null_space( vm[prefixvm( prefix,"constant-null-space" )].template as<bool>() ),
    M_showKSPMonitor( vm.count(prefixvm( prefix,"ksp-monitor" )) ),
    M_showSNESMonitor( vm.count(prefixvm( prefix,"snes-monitor" )) ),
    M_showKSPConvergedReason( vm.count(prefixvm( prefix,"ksp-converged-reason" )) ),
    M_showSNESConvergedReason( vm.count(prefixvm( prefix,"snes-converged-reason" )) )
{
}
template <typename T>
Backend<T>::~Backend()
{
    this->clear();
}
template <typename T>
void
Backend<T>::clear()
{
    if ( M_preconditioner )
        M_preconditioner->clear();
    LOG(INFO) << "Sending delete signal to all observers...\n";
    this->sendDeleteSignal();
    //this->clear ();
}
template <typename T>
typename Backend<T>::backend_ptrtype
Backend<T>::build( BackendType bt, WorldComm const& worldComm )
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
        return backend_ptrtype( new BackendPetsc<value_type>( worldComm ) );
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
template <typename T>
typename Backend<T>::backend_ptrtype
Backend<T>::build( po::variables_map const& vm, std::string const& prefix, WorldComm const& worldComm )
{

    std::string n = option( _name="backend" ).template as<std::string>();
    LOG(INFO) << "Loading backend " << n;
    BackendType bt;

    if ( ( n == "eigen_dense" ) || ( prefix.find("dd") != std::string::npos ) )
        bt = BACKEND_EIGEN_DENSE;

    else if ( n == "eigen" )
        bt = BACKEND_EIGEN;

    else if ( n == "petsc" )
        bt = BACKEND_PETSC;

    else if ( n == "trilinos" )
        bt = BACKEND_TRILINOS;

    else
    {

#if defined( FEELPP_HAS_PETSC_H )

        LOG(INFO) << "Falling back to backend petsc\n";
        bt = BACKEND_PETSC;
#else
        LOG(FATAL) << "Backend " << n << " not available";
#endif
    }

    // Build the appropriate solver
    switch ( bt )
    {
    case BACKEND_EIGEN:
    {
        return backend_ptrtype( new BackendEigen<value_type>( Environment::vm(), prefix, worldComm ) );
    }
    break;
    case BACKEND_EIGEN_DENSE:
    {
        return backend_ptrtype( new BackendEigen<value_type,1>( Environment::vm(), prefix, worldComm ) );
    }
    break;

#if defined ( FEELPP_HAS_PETSC_H )

    default:
    case BACKEND_PETSC:
    {
        return backend_ptrtype( new BackendPetsc<value_type>( Environment::vm(), prefix, worldComm ) );
    }
    break;
#endif
#if defined ( FEELPP_HAS_TRILINOS_EPETRA )

    case BACKEND_TRILINOS:
    {
#if defined ( FEELPP_HAS_TRILINOS_EPETRA )
        return backend_ptrtype( new BackendTrilinos( Environment::vm(), prefix, worldComm ) );
#else
        return backend_ptrtype();
#endif
    }
    break;
#endif
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

    MatrixStructure matStructInitial = this->precMatrixStructure();

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
        if (this->comm().globalRank() == this->comm().masterRank() )
            std::cout << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";
        LOG(INFO) << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";
        boost::tie( M_converged, M_iteration, M_residual ) = this->solve( A, P, x, b );

        //if ( !M_converged ) throw std::logic_error( "solver failed to converge" );
        if ( !M_converged ) std::cerr<< "Backend " << M_prefix << " : linear solver failed to converge" << std::endl;

        //stop();
    }

    this->setPrecMatrixStructure( matStructInitial );

    return boost::make_tuple( M_converged, M_iteration, M_residual );
}
template <typename T>
typename Backend<T>::nl_solve_return_type
Backend<T>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its,
                     bool reusePC, bool reuseJac )
{
    MatrixStructure matStructInitial = this->precMatrixStructure();

    M_nlsolver->setType( this->snesEnumType() );
    M_nlsolver->setPreconditionerType( this->pcEnumType() );
    M_nlsolver->setKspSolverType( this->kspEnumType() );
    M_nlsolver->setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );
    M_nlsolver->setShowSNESMonitor( this->showSNESMonitor() );
    M_nlsolver->setShowKSPMonitor( this->showKSPMonitor() );
    M_nlsolver->setShowKSPConvergedReason( this->showKSPConvergedReason() );
    M_nlsolver->setShowSNESConvergedReason( this->showSNESConvergedReason() );

    M_nlsolver->setNbItMax( this->maxIterationsSNES() );
    M_nlsolver->setRelativeResidualTol( this->rToleranceSNES() );
    M_nlsolver->setAbsoluteResidualTol( this->aToleranceSNES() );
    M_nlsolver->setAbsoluteSolutionTol( this->sToleranceSNES() );
    M_nlsolver->setRtoleranceKSP( this->rtoleranceKSPinSNES() );
    M_nlsolver->setAtoleranceKSP( this->aTolerance() );
    M_nlsolver->setDtoleranceKSP( this->dTolerance() );
    M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNES() );

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

        // if first time or rebuild prec at first newton step, need to get matStructInitial
        if ( reusePC && (!M_reusePrecIsBuild || M_reusePrecRebuildAtFirstNewtonStep) )
            {
                M_nlsolver->setPrecMatrixStructure( matStructInitial );
                M_reusePrecIsBuild=true;
            }
        else if ( reusePC ) M_nlsolver->setPrecMatrixStructure( SAME_PRECONDITIONER );

        // configure reusePC,reuseJac in non linear solver
        int typeReusePrec = 1,typeReuseJac = 1 ;
        if ( reusePC ) typeReusePrec = -1;
        if ( reuseJac ) typeReuseJac = -1;

        //M_nlsolver->setReuse( -2, -2 );
        //M_nlsolver->setReuse( -1, -2 );
        M_nlsolver->setReuse( typeReuseJac, typeReusePrec );

        //int maxIterationsReuseJac=10;
        M_nlsolver->setNbItMax( this->maxIterationsSNESReuse() );
        M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNESReuse() );

        // compute cst jacobian in case of quasi-newton!
        if ( reuseJac &&  (!M_reuseJacIsBuild || M_reuseJacRebuildAtFirstNewtonStep) ) { this->nlSolver()->jacobian( x, A );M_reuseJacIsBuild=true;}
    }

    auto ret = M_nlsolver->solve( A, x, b, tol, its );

    //std::cout << "[nlSolve] ret.first " << ret.first <<std::endl;
    if ( ret.first < 0 && ( reusePC || reuseJac ) )
    {
        if (this->comm().globalRank() == this->comm().masterRank() )
            std::cout << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";
        LOG(INFO) << "Backend "  << M_prefix << " reuse failed, rebuilding preconditioner...\n";

        // reset to initial solution
        x_save->close();
        *x=*x_save;

        //M_nlsolver->init();
        //M_nlsolver->setPreconditionerType( this->pcEnumType() );
        //M_nlsolver->setKspSolverType( this->kspEnumType() );
        M_nlsolver->setPrecMatrixStructure( matStructInitial/*SAME_NONZERO_PATTERN*/ );
        M_nlsolver->setNbItMax( this->maxIterationsSNES() );
        M_nlsolver->setMaxitKSP( this->maxIterationsKSPinSNES() );

        //M_nlsolver->setReuse( 1, -2 );
        M_nlsolver->setReuse( 1, 1 );

        this->nlSolver()->jacobian( x, A );

        // call solver which must execute with success
        auto ret2 = M_nlsolver->solve( A, x, b, tol, its );

        if ( ret2.first < 0 )
        {
            LOG(INFO) << "\n[backend] non-linear solver fail";
            //exit( 0 );
            std::cerr<< "Backend " << M_prefix << " : non-linear solver failed to converge" << std::endl;
        }

        return boost::make_tuple( ret2.first, its, tol );
    }

    this->setPrecMatrixStructure( matStructInitial );

    return boost::make_tuple( ret.first, its, tol );
}
template <typename T>
typename Backend<T>::nl_solve_return_type
Backend<T>::nlSolve( sparse_matrix_ptrtype& A,
                     vector_ptrtype& x,
                     vector_ptrtype& b,
                     const double tol, const int its )
{

    M_nlsolver->setType( this->snesEnumType() );
    M_nlsolver->setPreconditionerType( this->pcEnumType() );
    M_nlsolver->setKspSolverType( this->kspEnumType() );
    M_nlsolver->setMatSolverPackageType( this->matSolverPackageEnumType() );
    M_nlsolver->setPrecMatrixStructure( this->precMatrixStructure() );
    M_nlsolver->setShowSNESMonitor( this->showSNESMonitor() );
    M_nlsolver->setShowKSPMonitor( this->showKSPMonitor() );
    M_nlsolver->setShowKSPConvergedReason( this->showKSPConvergedReason() );
    M_nlsolver->setShowSNESConvergedReason( this->showSNESConvergedReason() );
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

    if ( ret.first < 0 )
    {
        LOG(ERROR) << "\n[backend] non-linear solver fail";
        LOG(ERROR) << "Backend " << M_prefix << " : non-linear solver failed to converge" << std::endl;
    }

    return boost::make_tuple( true, its, tol );
}
template <typename T>
typename Backend<T>::real_type
Backend<T>::dot( vector_type const& x, vector_type const& y ) const
{
    real_type localres = 0;

    for ( size_type i = 0; i < x.localSize(); ++i )
    {
        localres += x( i )*y( i );
    }

    real_type globalres=localres;
    mpi::all_reduce( M_worldComm.globalComm(), localres, globalres, std::plus<real_type>() );
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

template <typename T>
void
Backend<T>::reset()
{
    M_reusePC = false;
    M_totalSolveIter = 0.0;
    M_nUsePC = 0;
    //M_backend->set_maxiter( M_maxitKSP );

}

template<typename T>
SolverType
Backend<T>::kspEnumType() const
{
    if ( this->kspType()=="cg" )              return CG;

    else if ( this->kspType()=="cr" )         return CR;

    else if ( this->kspType()=="cgs" )        return CGS;

    else if ( this->kspType()=="bicg" )       return BICG;

    else if ( this->kspType()=="tcqmr" )      return TCQMR;

    else if ( this->kspType()=="tfqmr" )      return TFQMR;

    else if ( this->kspType()=="lsqr" )       return LSQR;

    else if ( this->kspType()=="bicgstab" )   return BICGSTAB;

    else if ( this->kspType()=="minres" )     return MINRES;

    else if ( this->kspType()=="gmres" )      return GMRES;

    else if ( this->kspType()=="richardson" ) return RICHARDSON;

    else if ( this->kspType()=="chebyshev" )  return CHEBYSHEV;

    else if ( this->kspType()=="preonly" )  return PREONLY;

    else return GMRES;

} // Backend::kspEnumType

template<typename T>
SolverNonLinearType
Backend<T>::snesEnumType() const
{
    if ( this->snesType() == "ls" || this->snesType() == "newtonls" )
        return SolverNonLinearType::LINE_SEARCH;
    else if ( this->snesType() == "tr" || this->snesType() == "newtontr" )
        return SolverNonLinearType::TRUST_REGION;
    else if ( this->snesType() == "nrichardson" )
        return SolverNonLinearType::NRICHARDSON;
    else if ( this->snesType() == "ksponly" )
        return SolverNonLinearType::NKSPONLY;
    else if ( this->snesType() == "vinewtonrsls" )
        return SolverNonLinearType::VINEWTONRSLS;
    else if ( this->snesType() == "vinewtonssls" )
        return SolverNonLinearType::VINEWTONRSTR;
    else if ( this->snesType() == "ngmres" )
        return SolverNonLinearType::NGMRES;
    else if ( this->snesType() == "qn" )
        return SolverNonLinearType::QN;
    else if ( this->snesType() == "shell" )
        return SolverNonLinearType::NSHELL;
    else if ( this->snesType() == "gs" )
        return SolverNonLinearType::GS;
    else if ( this->snesType() == "ncg" )
        return SolverNonLinearType::NCG;
    else if ( this->snesType() == "fas" )
        return SolverNonLinearType::FAS;
    else if ( this->snesType() == "ms" )
        return SolverNonLinearType::MS;
    else if ( this->snesType() == "nasm" )
        return SolverNonLinearType::NASM;
    else if ( this->snesType() == "anderson" )
        return SolverNonLinearType::ANDERSON;
    else if ( this->snesType() == "aspin" )
        return SolverNonLinearType::ASPIN;
    else
        return SolverNonLinearType::LINE_SEARCH;
}


template<typename T>
PreconditionerType
Backend<T>::pcEnumType() const
{

    if ( this->pcType()=="lu" )                return LU_PRECOND;

    else if ( this->pcType()=="ilu" )          return ILU_PRECOND;

    else if ( this->pcType()=="id" )           return IDENTITY_PRECOND;

    else if ( this->pcType()=="cholesky" )     return CHOLESKY_PRECOND;

    else if ( this->pcType()=="icc" )          return ICC_PRECOND;

    else if ( this->pcType()=="asm" )          return ASM_PRECOND;

    else if ( this->pcType()=="gasm" )          return GASM_PRECOND;

    else if ( this->pcType()=="jacobi" )       return JACOBI_PRECOND;

    else if ( this->pcType()=="block_jacobi" ) return BLOCK_JACOBI_PRECOND;

    else if ( this->pcType()=="bjacobi" )      return BLOCK_JACOBI_PRECOND;

    else if ( this->pcType()=="sor" )          return SOR_PRECOND;

    else if ( this->pcType()=="eisenstat" )    return EISENSTAT_PRECOND;

    else if ( this->pcType()=="shell" )        return SHELL_PRECOND;

    else if ( this->pcType()=="fieldsplit" )   return FIELDSPLIT_PRECOND;

    else if ( this->pcType()=="ml" )   return ML_PRECOND;

    else return LU_PRECOND;

} // Backend::pcEnumType

template<typename T>
FieldSplitType
Backend<T>::fieldSplitEnumType() const
{
    if ( this->fieldsplitType()=="additive" )  return ADDITIVE;

    else if ( this->fieldsplitType()=="multiplicative" )  return MULTIPLICATIVE;

    else if ( this->fieldsplitType()=="schur" )  return SCHUR;

    else return ADDITIVE;

} //Backend fieldSplitEnumType

template<typename T>
MatSolverPackageType
Backend<T>::matSolverPackageEnumType() const
{
    if ( this->pcFactorMatSolverPackageType()=="spooles" )           return MATSOLVER_SPOOLES;

    else if ( this->pcFactorMatSolverPackageType()=="superlu" )      return MATSOLVER_SUPERLU;

    else if ( this->pcFactorMatSolverPackageType()=="superlu-dist" ) return MATSOLVER_SUPERLU_DIST;

    else if ( this->pcFactorMatSolverPackageType()=="umfpack" )      return MATSOLVER_UMFPACK;

    else if ( this->pcFactorMatSolverPackageType()=="essl" )         return MATSOLVER_ESSL;

    else if ( this->pcFactorMatSolverPackageType()=="lusol" )        return MATSOLVER_LUSOL;

    else if ( this->pcFactorMatSolverPackageType()=="mumps" )        return MATSOLVER_MUMPS;

    else if ( this->pcFactorMatSolverPackageType()=="pastix" )       return MATSOLVER_PASTIX;

    else if ( this->pcFactorMatSolverPackageType()=="dscpack" )      return MATSOLVER_DSCPACK;

    else if ( this->pcFactorMatSolverPackageType()=="matlab" )       return MATSOLVER_MATLAB;

    else if ( this->pcFactorMatSolverPackageType()=="petsc" )        return MATSOLVER_PETSC;

    else if ( this->pcFactorMatSolverPackageType()=="plapack" )      return MATSOLVER_PLAPACK;

    else if ( this->pcFactorMatSolverPackageType()=="bas" )          return MATSOLVER_BAS;

    else return MATSOLVER_PETSC;
} // Backend matSolverPackageEnumType


/*
 * Explicit instantiations
 */
template class Backend<double>;


void updateBackendPreconditionerOptions( po::options_description & _options, std::string const& prefix )
{
    _options.add_options()
    ( prefixvm( prefix,"pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "lu" ), "type of preconditioners (lu, ilut, ilutp, diag, id,...)" )
    ( prefixvm( prefix,"sub-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "lu" ), "type of sub-preconditioners (lu, ilut, ilutp, diag, id,...)" )
    ( prefixvm( prefix,"pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "display preconditioner information" )
    ( prefixvm( prefix,"sub-pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "display sub-preconditioner information" )
    ( prefixvm( prefix,"constant-null-space" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "set the null space to be the constant values" )

    ( prefixvm( prefix,"pc-gasm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of gasm (basic, restrict, interpolate, none)" )
    ( prefixvm( prefix,"pc-gasm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of overlap levels" )
    ( prefixvm( prefix,"pc-asm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of asm (basic, restrict, interpolate, none)" )
    ( prefixvm( prefix,"pc-asm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of overlap levels" )
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
      "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
    ( prefixvm( prefix,"sub-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
      "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )

#else
    ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
      "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
    ( prefixvm( prefix,"sub-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
      "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif

    ( prefixvm( prefix,"ilu-threshold" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "threshold value for preconditioners" )
    ( prefixvm( prefix,"ilu-fillin" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "fill-in level value for preconditioners" )
    ( prefixvm( prefix,"pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu" )
    ( prefixvm( prefix,"sub-pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu as a sub-preconditioner" )
    ( prefixvm( prefix,"pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ), "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )
    ( prefixvm( prefix,"sub-pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ), "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )
        ;
}

/**
 * \return the command lines options of the petsc backend
 */
po::options_description backend_options( std::string const& prefix )
{
    po::options_description _options( "Linear and NonLinear Solvers Backend " + prefix + " options" );
    _options.add_options()
        // solver options
        ( prefixvm( prefix,"backend" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ), "backend type: petsc, eigen, eigen_dense" )
        ( prefixvm( prefix,"backend.verbose" ).c_str(), Feel::po::value<bool>()->default_value( false ), "set the backend to be verbose" )
        ( prefixvm( prefix,"ksp-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-13 ), "relative tolerance" )
        ( prefixvm( prefix,"ksp-atol" ).c_str(), Feel::po::value<double>()->default_value( 1e-50 ), "absolute tolerance" )
        ( prefixvm( prefix,"ksp-dtol" ).c_str(), Feel::po::value<double>()->default_value( 1e5 ), "divergence tolerance" )
        ( prefixvm( prefix,"ksp-maxit" ).c_str(), Feel::po::value<size_type>()->default_value( 1000 ), "maximum number of iterations" )
        ( prefixvm( prefix,"ksp-maxit-reuse" ).c_str(), Feel::po::value<size_type>(), "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,"reuse-jac" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reuse jacobian" )
        ( prefixvm( prefix,"reuse-jac.rebuild-at-first-newton-step" ).c_str(), Feel::po::value<bool>()->default_value( true ), "rebuild jacobian at each Newton when reuse jacobian" )
        ( prefixvm( prefix,"reuse-prec" ).c_str(), Feel::po::value<bool>()->default_value( false ), "reuse preconditioner" )
        ( prefixvm( prefix,"reuse-prec.rebuild-at-first-newton-step" ).c_str(), Feel::po::value<bool>()->default_value( false ), "rebuild preconditioner at each Newton when reuseprec" )

        ( prefixvm( prefix,"export-matlab" ).c_str(), Feel::po::value<std::string>()->default_value( "" ), "export matrix/vector to matlab, default empty string means no export, other string is used as prefix" )

        ( prefixvm( prefix,"ksp-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Prints the KSP data structure" )
        ( prefixvm( prefix,"ksp-type" ).c_str(), Feel::po::value<std::string>()->default_value( "gmres" ), "cg, bicgstab, gmres" )
        ( prefixvm( prefix,"ksp-monitor" ).c_str(), Feel::po::value<bool>()->default_value( false ) , "monitor ksp" )
        ( prefixvm( prefix,"ksp-converged-reason" ).c_str() , "converged reason ksp" )

        ( prefixvm( prefix,"snes-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Prints the SNES data structure" )
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
        ( prefixvm( prefix,"snes-type" ).c_str(), Feel::po::value<std::string>()->default_value( SNESNEWTONLS ), "Set the SNES solver" )
#else
        ( prefixvm( prefix,"snes-type" ).c_str(), Feel::po::value<std::string>()->default_value( SNESLS ), "Set the SNES solver" )
#endif
        ( prefixvm( prefix,"snes-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-8 ), "relative tolerance" )
        ( prefixvm( prefix,"snes-atol" ).c_str(), Feel::po::value<double>()->default_value( 1e-50 ), "absolute tolerance" )
        ( prefixvm( prefix,"snes-stol" ).c_str(), Feel::po::value<double>()->default_value( 1e-8 ), "step length tolerance" )
        ( prefixvm( prefix,"snes-maxit" ).c_str(), Feel::po::value<size_type>()->default_value( 50 ), "maximum number of iterations" )
        ( prefixvm( prefix,"snes-maxit-reuse" ).c_str(), Feel::po::value<size_type>(), "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,"snes-ksp-maxit" ).c_str(), Feel::po::value<size_type>()->default_value( 1000 ), "maximum number of iterations" )
        ( prefixvm( prefix,"snes-ksp-maxit-reuse" ).c_str(), Feel::po::value<size_type>(), "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,"snes-ksp-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-5 ), "relative tolerance" )
        ( prefixvm( prefix,"snes-monitor" ).c_str(), Feel::po::value<bool>()->default_value( false ) , "monitor snes" )
        ( prefixvm( prefix,"snes-converged-reason" ).c_str() , "converged reason snes" )


        // preconditioner options
        ( prefixvm( prefix,"pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Prints the PC data structure" )
        ( prefixvm( prefix,"pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "lu" ), "type of preconditioners (lu, ilut, ilutp, diag, id,...)" )
        ( prefixvm( prefix,"sub-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "lu" ), "type of sub-preconditioners (lu, ilut, ilutp, diag, id,...)" )
        ( prefixvm( prefix,"sub-pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "display sub-preconditioner information" )
        ( prefixvm( prefix,"constant-null-space" ).c_str(), Feel::po::value<bool>()->default_value( 0 ), "set the null space to be the constant values" )

        ( prefixvm( prefix,"pc-gasm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of gasm (basic, restrict, interpolate, none)" )
        ( prefixvm( prefix,"pc-gasm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "number of overlap levels" )
        ( prefixvm( prefix,"pc-asm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of asm (basic, restrict, interpolate, none)" )
        ( prefixvm( prefix,"pc-asm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "number of overlap levels" )
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
        ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
        ( prefixvm( prefix,"sub-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )

#else
        ( prefixvm( prefix,"pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
        ( prefixvm( prefix,"sub-pc-factor-mat-solver-package-type" ).c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif

        ( prefixvm( prefix,"ilu-threshold" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "threshold value for preconditioners" )
        ( prefixvm( prefix,"ilu-fillin" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "fill-in level value for preconditioners" )
        ( prefixvm( prefix,"pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu" )
        ( prefixvm( prefix,"sub-pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu as a sub-preconditioner" )
        ( prefixvm( prefix,"pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ), "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )
        ( prefixvm( prefix,"sub-pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ), "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )

        // solver control options
        ( prefixvm( prefix,"gmres-restart" ).c_str(), Feel::po::value<int>()->default_value( 30 ), "number of iterations before solver restarts (gmres)" )
        ( prefixvm( prefix,"ksp-verbose" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "(=0,1,2) print solver iterations" )
        ( prefixvm( prefix,"fieldsplit-type" ).c_str(), Feel::po::value<std::string>()->default_value( "additive" ), "type of fieldsplit (additive, multiplicative, schur)" )
        ( prefixvm( prefix,"fieldsplit-schur-fact-type" ).c_str(), Feel::po::value<std::string>()->default_value( "full" ), "type of schur factorization (diag, lower, upper, full)" )
        ;

    for ( uint16_type i=0; i<5; ++i )
    {
        std::string prefixfieldsplit = ( boost::format( "%1%fieldsplit-%2%" ) %prefixvm( prefix,"" ) %i ).str();

        _options.add_options()
            ( prefixvm( prefixfieldsplit,"pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( (i==0)?"lu":"none" ), "type of fieldsplit preconditioners" )
            ( prefixvm( prefixfieldsplit,"sub-pc-type" ).c_str(), Feel::po::value<std::string>()->default_value( "lu" ), "type of fieldsplit preconditioners" )
            ( prefixvm( prefixfieldsplit,"ksp-type").c_str(), Feel::po::value<std::string>()->default_value( /*"gmres"*/"preonly" ), "type of fieldsplit solver" )
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
            (  prefixvm( prefixfieldsplit,"pc-factor-mat-solver-package-type").c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
               "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
            (  prefixvm( prefixfieldsplit,"sub-pc-factor-mat-solver-package-type").c_str(), Feel::po::value<std::string>()->default_value( "mumps" ),
               "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#else
            ( prefixvm( prefixfieldsplit,"pc-factor-mat-solver-package-type").c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
              "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
            ( prefixvm( prefixfieldsplit,"sub-pc-factor-mat-solver-package-type").c_str(), Feel::po::value<std::string>()->default_value( "petsc" ),
              "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif

            ( prefixvm( prefixfieldsplit,"pc-gasm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of gasm (basic, restrict, interpolate, none)" )
            ( prefixvm( prefixfieldsplit,"pc-gasm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "number of overlap levels" )
            ( prefixvm( prefixfieldsplit,"pc-asm-type" ).c_str(), Feel::po::value<std::string>()->default_value( "restrict" ), "type of asm (basic, restrict, interpolate, none)" )
            ( prefixvm( prefixfieldsplit,"pc-asm-overlap" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "number of overlap levels" )

            ( prefixvm( prefixfieldsplit,"pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu as a sub-preconditioner" )
            ( prefixvm( prefixfieldsplit,"sub-pc-factor-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "Sets the number of levels of fill to use for ilu as a sub-preconditioner" )
            ( prefixvm( prefixfieldsplit,"pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ),
              "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )
            ( prefixvm( prefixfieldsplit,"sub-pc-factor-fill" ).c_str(), Feel::po::value<double>()->default_value( 6 ),
              "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )

            ( prefixvm( prefixfieldsplit,"pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "display sub-preconditioner information" )
            ( prefixvm( prefixfieldsplit,"sub-pc-view" ).c_str(), Feel::po::value<bool>()->default_value( false ), "display sub-preconditioner information" )
            ;
    }

    // mg or ml preconditioner
    _options.add_options()
        ( prefixvm( prefix,"pc-mg-levels" ).c_str(), Feel::po::value<int>()->default_value( 3 ), "number of levels including finest" )
        ( prefixvm( prefix,"pc-mg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "kaskade" ), "Determines the form of multigrid to use: multiplicative, additive, full, kaskade " )
        ( prefixvm( prefix,"pc-mg-smoothdown" ).c_str(), Feel::po::value<int>()->default_value( 1 ), "number of smoothing steps before applying restriction operator" )
        ( prefixvm( prefix,"pc-ml-reuse-interpolation" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Reuse the interpolation operators when possible (cheaper, weaker when matrix entries change a lot)" )
        ( prefixvm( prefix,"pc-ml-keep-agg-info" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Allows the preconditioner to be reused, or auxilliary matrices to be generated" )
        ( prefixvm( prefix,"pc-ml-reusable" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Store intermedaiate data structures so that the multilevel hierarchy is reusable" )
        ( prefixvm( prefix,"pc-ml-old-hierarchy" ).c_str(), Feel::po::value<bool>()->default_value( false ), "Use old routine to generate hierarchy" )
        ;

    std::string prefixMGCoarse = ( boost::format( "%1%mg-coarse" ) %prefixvm( prefix,"" ) ).str();
    updateBackendPreconditionerOptions( _options, prefixMGCoarse );

    for ( uint16_type i=1; i<6; ++i )
    {
        std::string prefixMGLevels = ( boost::format( "%1%mg-levels%2%" ) %prefixvm( prefix,"" ) %i ).str();
        updateBackendPreconditionerOptions( _options, prefixMGLevels );

        _options.add_options()
            ( prefixvm( prefixMGLevels,"ksp-type" ).c_str(), Feel::po::value<std::string>()->default_value( "gmres" ), "cg, bicgstab, gmres" )
            ( prefixvm( prefixMGLevels,"ksp-monitor" ).c_str(), Feel::po::value<bool>()->default_value( false ) , "monitor ksp" )
            ( prefixvm( prefixMGLevels,"ksp-rtol" ).c_str(), Feel::po::value<double>()->default_value( 1e-3 ), "relative tolerance" )
            ( prefixvm( prefixMGLevels,"ksp-maxit" ).c_str(), Feel::po::value<size_type>()->default_value( 50 ), "maximum number of iterations" )
            ;
    }

    return _options;
}



} // namespace Feel

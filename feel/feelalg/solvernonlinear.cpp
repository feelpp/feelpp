/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-04-17

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file solvernonlinear.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-04-17
 */
#define FEELPP_INSTANTIATE_SOLVERNONLINEAR

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/solvernonlinearpetsc.hpp>
#include <feel/feelalg/solvernonlineartrilinos.hpp>

namespace Feel
{
template <typename T>
inline
SolverNonLinear<T>::SolverNonLinear ( std::string const& prefix, WorldComm const& worldComm )
    :
    residual        ( 0 ),
    jacobian        ( 0 ),
    matvec          ( 0 ),
    M_prefix( prefix ),
    M_worldComm( worldComm ),
    M_is_initialized ( false ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_snl_type( snesTypeConvertStrToEnum( soption(_prefix=prefix,_name="snes-type" ) ) ),
    M_kspSolver_type( GMRES ),
    M_preconditioner_type( LU_PRECOND ),
    M_preconditioner(),
    M_matSolverPackage_type( MATSOLVER_PETSC ),
    M_relativeResidualTol( 0 ),
    M_absoluteResidualTol( 0 ),
    M_absoluteSolutionTol( 0 ),
    M_nbItMax( 0 ),
    M_reuse_jac( 0 ),
    M_reuse_prec( 0 ),
    M_showKSPMonitor( boption(_prefix=prefix,_name="ksp-monitor") ),
    M_showSNESMonitor( boption(_prefix=prefix,_name="snes-monitor" ) ),
    M_showKSPConvergedReason( Environment::vm().count(prefixvm( prefix,"ksp-converged-reason" )) ),
    M_showSNESConvergedReason( Environment::vm().count(prefixvm( prefix,"snes-converged-reason" )) ),
    M_viewSNESInfo( boption( _prefix=prefix, _name="snes-view" ) ),
    M_rtoleranceKSP( 1e-13 ),
    M_dtoleranceKSP( 1e5 ),
    M_atoleranceKSP( 1e-50 ),
    M_maxitKSP( 1000 )
{
}

template <typename T>
inline
SolverNonLinear<T>::SolverNonLinear ( SolverNonLinear const& snl )
    :
    residual        ( snl.residual ),
    jacobian        ( snl.jacobian ),
    matvec          ( snl.matvec ),
    M_prefix( snl.M_prefix ),
    M_worldComm( snl.M_worldComm ),
    M_is_initialized( snl.M_is_initialized ),
    M_prec_matrix_structure( snl.M_prec_matrix_structure ),
    M_snl_type( snl.M_snl_type ),
    M_kspSolver_type( snl.M_kspSolver_type ),
    M_preconditioner_type( snl.M_preconditioner_type ),
    M_preconditioner( snl.M_preconditioner ),
    M_matSolverPackage_type( snl.M_matSolverPackage_type ),
    M_relativeResidualTol( snl.M_relativeResidualTol ),
    M_absoluteResidualTol( snl.M_absoluteResidualTol ),
    M_absoluteSolutionTol( snl.M_absoluteSolutionTol ),
    M_nbItMax( snl.M_nbItMax ),
    M_reuse_jac( snl.M_reuse_jac ),
    M_reuse_prec( snl.M_reuse_prec ),
    M_showKSPMonitor( snl.M_showKSPMonitor ),
    M_showSNESMonitor( snl.M_showSNESMonitor ),
    M_showKSPConvergedReason( snl.M_showKSPConvergedReason ), M_showSNESConvergedReason( snl.M_showSNESConvergedReason ),
    M_viewSNESInfo( snl.M_viewSNESInfo ),
    M_rtoleranceKSP( snl.M_rtoleranceKSP ),
    M_dtoleranceKSP( snl.M_dtoleranceKSP ),
    M_atoleranceKSP( snl.M_atoleranceKSP ),
    M_maxitKSP( snl.M_maxitKSP )
{
}



template <typename T>
inline
SolverNonLinear<T>::~SolverNonLinear ()
{
    this->clear ();
}

template <typename T>
boost::shared_ptr<SolverNonLinear<T> >
SolverNonLinear<T>::build( po::variables_map const& vm, std::string const& prefix, WorldComm const& worldComm )
{
    return build( prefix,worldComm );
}
template <typename T>
boost::shared_ptr<SolverNonLinear<T> >
SolverNonLinear<T>::build( std::string const& prefix, WorldComm const& worldComm )
{
    return build( soption(_name="backend"),prefix,worldComm );
}
template <typename T>
boost::shared_ptr<SolverNonLinear<T> >
SolverNonLinear<T>::build( std::string const& kind, std::string const& prefix, WorldComm const& worldComm )
{
    SolverPackage solver_package=SOLVER_INVALID_PACKAGE;

    if ( kind == "petsc" )
    {
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_PETSC;
#endif
    }

    else if ( kind == "trilinos" )
    {
#if defined( FEELPP_HAS_TRILINOS )
        solver_package = SOLVERS_TRILINOS;
#endif
    }

    else
    {
        LOG(INFO) << "[SolverNonLinear] solver " << kind << " not available\n";
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_PETSC;
#endif
    }

    // Build the appropriate solver
    switch ( solver_package )
    {
#if defined( FEELPP_HAS_PETSC )

    case SOLVERS_PETSC:
    {


        solvernonlinear_ptrtype ap( new SolverNonLinearPetsc<T>( prefix,worldComm ) );
        return ap;
    }
    break;
#endif

#if defined( FEELPP_HAS_TRILINOS )

    case SOLVERS_TRILINOS:
    {
        solvernonlinear_ptrtype ap( new SolverNonLinearTrilinos<T>( prefix,worldComm )  );
        return ap;
    }
    break;
#endif

    default:
        std::cerr << "ERROR:  Unrecognized NonLinear solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

    return solvernonlinear_ptrtype();
}
template <>
boost::shared_ptr<SolverNonLinear<std::complex<double>> >
SolverNonLinear<std::complex<double>>::build( std::string const& kind, std::string const& prefix, WorldComm const& worldComm )
{
    SolverPackage solver_package=SOLVER_INVALID_PACKAGE;

    if ( kind == "petsc" )
    {
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_PETSC;
#endif
    }

    else if ( kind == "trilinos" )
    {
#if defined( FEELPP_HAS_TRILINOS )
        solver_package = SOLVERS_TRILINOS;
#endif
    }

    else
    {
        LOG(INFO) << "[SolverNonLinear] solver " << kind << " not available\n";
        LOG(INFO) << "[Backend] use fallback  gmm\n";
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_PETSC;
#endif
    }

    // Build the appropriate solver
    switch ( solver_package )
    {
#if defined( FEELPP_HAS_PETSC ) && defined (PETSC_HAS_COMPLEX_SUPPORT)

    case SOLVERS_PETSC:
    {


        solvernonlinear_ptrtype ap( new SolverNonLinearPetsc<T>( prefix,worldComm ) );
        return ap;
    }
    break;
#endif

#if defined( FEELPP_HAS_TRILINOS )

    case SOLVERS_TRILINOS:
    {
        solvernonlinear_ptrtype ap( new SolverNonLinearTrilinos<T>( prefix,worldComm )  );
        return ap;
    }
    break;
#endif

    default:
        std::cerr << "ERROR:  Unrecognized NonLinear solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

    return solvernonlinear_ptrtype();
}

template <typename T>
boost::shared_ptr<SolverNonLinear<T> >
SolverNonLinear<T>::build( SolverPackage solver_package, WorldComm const& worldComm )
{
#if defined( FEELPP_HAS_PETSC )

    if ( solver_package != SOLVERS_PETSC )
    {
        LOG(INFO) << "[SolverNonLinear] solver " << solver_package << " not available\n";
        LOG(INFO) << "[Backend] use fallback  petsc: " << SOLVERS_PETSC << "\n";
        solver_package = SOLVERS_PETSC;
    }

    // Build the appropriate solver
    switch ( solver_package )
    {

    case SOLVERS_PETSC:
    {

#if defined( FEELPP_HAS_PETSC )
        solvernonlinear_ptrtype ap(new SolverNonLinearPetsc<T>( "",worldComm ) );
        return ap;
#else
        std::cerr << "PETSc is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver PETSc package" );
#endif
    }
    break;

    case SOLVERS_TRILINOS:
    {
#if defined( FEELPP_HAS_TRILINOS )
        solvernonlinear_ptrtype ap( new SolverNonLinearTrilinos<T>( "",worldComm ) );
        return ap;
#else
        std::cerr << "Trilinos NOX is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver NOX package" );
#endif
    }
    break;

    default:
        std::cerr << "ERROR:  Unrecognized NonLinear solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

#endif
    return solvernonlinear_ptrtype();
}

template <>
boost::shared_ptr<SolverNonLinear<std::complex<double>> >
SolverNonLinear<std::complex<double>>::build( SolverPackage solver_package, WorldComm const& worldComm )
{
#if defined( FEELPP_HAS_PETSC )

    if ( solver_package != SOLVERS_PETSC )
    {
        LOG(INFO) << "[SolverNonLinear] solver " << solver_package << " not available\n";
        LOG(INFO) << "[Backend] use fallback  petsc: " << SOLVERS_PETSC << "\n";
        solver_package = SOLVERS_PETSC;
    }

    // Build the appropriate solver
    switch ( solver_package )
    {

    case SOLVERS_PETSC:
    {

#if defined( FEELPP_HAS_PETSC ) && defined (PETSC_HAS_COMPLEX_SUPPORT )
        solvernonlinear_ptrtype ap(new SolverNonLinearPetsc<T>( "",worldComm ) );
        return ap;
#else
        std::cerr << "PETSc is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver PETSc package" );
#endif
    }
    break;

    case SOLVERS_TRILINOS:
    {
#if defined( FEELPP_HAS_TRILINOS )
        solvernonlinear_ptrtype ap( new SolverNonLinearTrilinos<T>( "",worldComm ) );
        return ap;
#else
        std::cerr << "Trilinos NOX is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver NOX package" );
#endif
    }
    break;

    default:
        std::cerr << "ERROR:  Unrecognized NonLinear solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

#endif
    return solvernonlinear_ptrtype();
}


/*
 * Explicit instantiations
 */
template class SolverNonLinear<double>;
template class SolverNonLinear<std::complex<double>>;

/**
 * \return the command lines options of the petsc backend
 */
po::options_description nlsolver_options()
{
    po::options_description _options( "Non Linear Solver options" );
    //_options.add_options()
    // solver options
    //("backend", Feel::po::value<std::string>()->default_value( "petsc" ), "nonlinear solver type: petsc")
    //("backend", Feel::po::value<std::string>()->default_value( "trilinos" ), "nonlinear solver type: trilinos NOX")
    //;
    return _options;
}

}

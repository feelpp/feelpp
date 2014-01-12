/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-04

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file solvereigen.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-04
 */
#include <boost/algorithm/string.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelalg/solvereigenslepc.hpp>

namespace Feel
{
template <typename T>
SolverEigen<T>::SolverEigen()
    :
    M_eigen_solver_type    ( ARNOLDI ),
    M_eigen_problem_type   ( NHEP ),
    M_position_of_spectrum ( LARGEST_MAGNITUDE ),
    M_spectral_transform   ( SHIFT ),
    M_is_initialized       ( false ),
    M_nev( 1 ),
    M_ncv( 3 )
{
}

template <typename T>
SolverEigen<T>::SolverEigen( po::variables_map const& vm, std::string const& prefix )
    :
    M_prefix( ( !prefix.empty() && !boost::algorithm::ends_with( prefix, "-" ) )? ( prefix+"-" ):prefix ),
    M_eigen_solver_type    ( ( EigenSolverType )vm[M_prefix+"solvereigen.solver-type"].template as<int>() ),
    M_eigen_problem_type   ( ( EigenProblemType )vm[M_prefix+"solvereigen.problem-type"].template as<int>() ),
    M_position_of_spectrum ( ( PositionOfSpectrum )vm[M_prefix+"solvereigen.position"].template as<int>() ),
    M_spectral_transform   ( SHIFT ),
    M_is_initialized       ( false ),
    M_nev( vm[M_prefix+"solvereigen.nev"].template as<int>() ),
    M_ncv( vm[M_prefix+"solvereigen.ncv"].template as<int>() )
{
}

template <typename T>
SolverEigen<T>::SolverEigen( SolverEigen const& eis )
    :
    M_prefix( eis.M_prefix ),
    M_eigen_solver_type    ( eis.M_eigen_solver_type ),
    M_eigen_problem_type   ( eis.M_eigen_problem_type ),
    M_position_of_spectrum ( eis.M_position_of_spectrum ),
    M_spectral_transform   ( eis.M_spectral_transform ),
    M_is_initialized       ( eis.M_is_initialized ),
    M_nev( eis.M_nev ),
    M_ncv( eis.M_ncv ),
    M_maxit( eis.M_maxit ),
    M_tolerance( eis.M_tolerance ),
    M_mapRow( eis.M_mapRow ),
    M_mapCol( eis.M_mapCol )
{
}

template <typename T>
SolverEigen<T>::~SolverEigen()
{
    this->clear ();
}
template <typename T>
boost::shared_ptr<SolverEigen<T> >
SolverEigen<T>::build( const SolverPackage solver_package )
{
    // Build the appropriate solver
    switch ( solver_package )
    {

    case SOLVERS_SLEPC:
    {

#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
        solvereigen_ptrtype ap( new SolverEigenSlepc<T> );
        return ap;
#else
        std::cerr << "Slepc is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver slepc package" );
#endif
    }
    break;

    default:
        std::cerr << "ERROR:  Unrecognized eigen solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

    return solvereigen_ptrtype();
}

template <typename T>
boost::shared_ptr<SolverEigen<T> >
SolverEigen<T>::build( po::variables_map const& vm, std::string const& prefix )
{
    SolverPackage solver_package = SOLVERS_SLEPC;

    if ( vm["backend"].template as<std::string>() == "petsc" )
    {
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_SLEPC;
#endif
    }

    else if ( vm["backend"].template as<std::string>() == "trilinos" )
    {
#if defined( FEELPP_HAS_TRILINOS )
        solver_package = SOLVERS_TRILINOS;
#endif
    }

    else
    {
        LOG(INFO) << "[SolverNonLinear] solver " << vm["backend"].template as<std::string>() << " not available\n";
        LOG(INFO) << "[Backend] use fallback  gmm\n";
#if defined( FEELPP_HAS_PETSC )
        solver_package = SOLVERS_PETSC;
#endif
    }

    // Build the appropriate solver
    switch ( solver_package )
    {

    case SOLVERS_SLEPC:
    {

#if defined( FEELPP_HAS_SLEPC ) && defined( FEELPP_HAS_PETSC_H )
        solvereigen_ptrtype ap( new SolverEigenSlepc<T>( vm, prefix ) );
        return ap;
#else
        std::cerr << "Slepc is not available/installed" << std::endl;
        throw std::invalid_argument( "invalid solver slepc package" );
#endif
    }
    break;

    default:
        std::cerr << "ERROR:  Unrecognized eigen solver package: "
                  << solver_package
                  << std::endl;
        throw std::invalid_argument( "invalid solver package" );
    }

    return solvereigen_ptrtype();
}


/*
 * Explicit instantiations
 */
template class SolverEigen<double>;
}

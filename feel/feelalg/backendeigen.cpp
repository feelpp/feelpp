/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-05-22

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file backendeigen.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-05-22
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backendeigen.hpp>

namespace Feel
{
// -- CONSTRUCTOR --
template<typename T>
BackendEigen<T>::BackendEigen( WorldComm const& )
    :
    super()
{}

template<typename T>
BackendEigen<T>::BackendEigen( po::variables_map const& vm, std::string const& prefix, WorldComm const&  )
    :
    super( vm, prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";
}


template<typename T>
typename BackendEigen<T>::solve_return_type
BackendEigen<T>::solve( sparse_matrix_type const& _A,
                      vector_type& _x,
                      const vector_type& _b )
{
    bool reusePC = ( this->precMatrixStructure() == SAME_PRECONDITIONER );

    eigen_sparse_matrix_type const& A( dynamic_cast<eigen_sparse_matrix_type const&>( _A ) );
    eigen_vector_type      & x( dynamic_cast<eigen_vector_type      &>( _x ) );
    eigen_vector_type const& b( dynamic_cast<eigen_vector_type const&>( _b ) );

} // BackendEigen::solve

//
// Instantiation
//
template class BackendEigen<double>;

/**
 * \return the command lines options of the eigen backend
 */
po::options_description backendeigen_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";


    po::options_description _options( "BackendEigen " + prefix + " solver options" );
#if 0
    _options.add_options()
    // solver options
    ( ( _prefix+"eigen-solver-type" ).c_str(),
      Feel::po::value<std::string>()->default_value( defaults.solver_type ),
      "umfpack, superlu, cg, bicgstab, gmres" )

    // preconditioner options
    ( ( _prefix+"eigen-pc-type" ).c_str(),
      Feel::po::value<std::string>()->default_value( defaults.pc_type ),
      "ilut, ilutp, diag, id" )
    ( ( _prefix+"eigen-threshold" ).c_str(),
      Feel::po::value<double>()->default_value( defaults.threshold ),
      "threshold value for preconditioners" )
    ( ( _prefix+"eigen-fillin" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.fillin ),
      "fill-in level value for preconditioners" )

    // solver control options
    ( ( _prefix+"eigen-restart" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.restart ),
      "number of iterations before solver restarts (gmres)" )
    ( ( _prefix+"eigen-verbose" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.verbose ),
      "(=0,1,2) print solver iterations" )
    ( ( _prefix+"eigen-maxiter" ).c_str(),
      Feel::po::value<int>()->default_value( defaults.maxiter ),
      "set maximum number of iterations" )
    ( ( _prefix+"eigen-tolerance" ).c_str(),
      Feel::po::value<double>()->default_value( defaults.tolerance ),
      "set solver tolerance" )
    ;
#endif
    return _options;
}

} // Feel

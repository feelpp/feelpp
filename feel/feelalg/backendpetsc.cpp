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
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backendpetsc.hpp>

namespace Feel
{

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



} // Feel

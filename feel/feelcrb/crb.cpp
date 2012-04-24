/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-11-24

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file opuscrb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-11-24
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbscm.hpp>

namespace Feel
{
Feel::po::options_description
crbOptions( std::string const& prefix )
{
    Feel::po::options_description crboptions( "CRB Options" );
    crboptions.add_options()
    ( "crb.output-index"   , Feel::po::value<int>()->default_value( 0 ), "output index (default is right hand side = 0)" )
    ( "crb.sampling-size"   , Feel::po::value<int>()->default_value( 1000 ), "Offline  sampling size " )
    ( "crb.error-max"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
    ( "crb.online-tolerance"   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
    ( "crb.dimension-max"   , Feel::po::value<int>()->default_value( 50 ),       "Offline  max WN size" )
    ( "crb.error-type"   , Feel::po::value<int>()->default_value( ( int )CRB_RESIDUAL ),       "CRB error type to be computed" )
    ( "crb.factor"   , Feel::po::value<int>()->default_value( -1 ),  "factor useful to estimate error by empirical method" )
    ( "crb.Nm"   , Feel::po::value<int>()->default_value( 1 ),       "Offline  number of modes per mu (for the POD) " )
    ( "crb.check.rb"   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
    ( "crb.check.gs"   , Feel::po::value<int>()->default_value( 0 ),       "check Gram-Schmidt orthonormalisation" )
    ( "crb.check.residual"   , Feel::po::value<int>()->default_value( 0 ),  "check residual" )
    ( "crb.reuse-prec"   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
    ( "crb.orthonormalize-primal" , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
    ( "crb.orthonormalize-dual" , Feel::po::value<bool>()->default_value( 1 ), "orthonormalize or not " )
    ( "crb.solve-dual-problem" , Feel::po::value<bool>()->default_value( 1 ), "solve or not the dual problem (this bool will be ignored if error-type=CRB_RESIDUAL) " )
    ( "crb.visualize-basis" , Feel::po::value<bool>()->default_value( 0 ), "visualize elements of the reduced basis " )
    ( "crb.save-output-behavior" , Feel::po::value<bool>()->default_value( 0 ), "save output behavior in time" )
    ( "crb.seek-mu-in-complement" , Feel::po::value<bool>()->default_value( 1 ), "during the offline basis construction, see mu in M the complement of Wn" )
    ( "crb.rebuild-database" , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
    ( "crb.show-mu-selection" , Feel::po::value<bool>()->default_value( 0 ), " show mu selection during offline step to build RB space" )
    ( "crb.compute-variance" , Feel::po::value<bool>()->default_value( 0 ), " if true the output is the variance and not l(v)" )
    ;



    crboptions
    .add( crbSCMOptions() );

    return crboptions;
}
}


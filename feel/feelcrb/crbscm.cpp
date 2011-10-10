/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-08-10

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
   \file crbscm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-08-10
 */

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/crbscm.hpp>

namespace Feel
{
Feel::po::options_description
crbSCMOptions( std::string const& prefix )
{
    Feel::po::options_description crbscmoptions("CRB SCM Options");
    crbscmoptions.add_options()
        ("crb.scm.sampling-size"   , Feel::po::value<int>()->default_value( 100 ),       "Offline SCM sampling size " )
        ("crb.scm.tol"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline SCM tolerance" )
        ("crb.scm.iter-max"   , Feel::po::value<int>()->default_value( 10 ),       "Offline SCM max number of K iterations" )
        ("crb.scm.Mplus" , Feel::po::value<int>()->default_value( 10 ),       "M+ value" )
        ("crb.scm.Malpha" , Feel::po::value<int>()->default_value( 10 ),       "M_alpha value" )
        ("crb.scm.level" , Feel::po::value<int>()->default_value( 1 ),       "level for recursion in lower bound computations" )
        ("crb.scm.strategy" , Feel::po::value<int>()->default_value( 2 ),       "scm strategy (0=patera, 1=maday, 2=prudhomme" )
        ;



    crbscmoptions
        .add( solvereigen_options( "crb.scm" ) );

    return crbscmoptions;
  }
}

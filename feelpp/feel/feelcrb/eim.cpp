/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-30

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file eim.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-05-30
 */

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/eim.hpp>


namespace Feel
{
Feel::po::options_description
eimOptions( std::string const& prefix )
{
    Feel::po::options_description eimoptions( "EIM Options" );
    eimoptions.add_options()
        ( "eim.sampling-size"   , Feel::po::value<int>()->default_value( 30 ), "Offline  sampling size " )
        ( "eim.error-max"   , Feel::po::value<double>()->default_value( 1e-6 ),       "Offline  tolerance" )
        ( "eim.online-tolerance"   , Feel::po::value<double>()->default_value( 1e-2 ),       "Online  tolerance" )
        ( "eim.dimension-max"   , Feel::po::value<int>()->default_value( 50 ),       "Offline  max WN size" )
        //( "eim.error-type"   , Feel::po::value<int>()->default_value( ( int )EIM_RESIDUAL ),       "EIM error type to be computed" )
        ( "eim.check.rb"   , Feel::po::value<int>()->default_value( 0 ),       "check reduced basis" )
        ( "eim.norm-used-for-residual" , Feel::po::value<std::string>()->default_value( "Linfty" ), "name of the norm used : Linfty - L2 - LinftyVec" )
        ( "eim.check.residual"   , Feel::po::value<int>()->default_value( 0 ),  "check residual" )
        ( "eim.reuse-prec"   , Feel::po::value<bool>()->default_value( 0 ),       "reuse or not the preconditioner" )
        ( "eim.rebuild-database" , Feel::po::value<bool>()->default_value( 0 ), "rebuild database (if it already exists)" )
        ( "eim.cvg-study" , Feel::po::value<bool>()->default_value( 0 ), "for convergence study" )
        ( "eim.compute-error-with-truth-expression" , Feel::po::value<bool>()->default_value( true ), "compute the error with the truth expression ( not its projection ) if true" )
        ( "eim.use-dimension-max-functions" , Feel::po::value<bool>()->default_value( 0 ), "force to use dimension-max basis functions" )
        ( "eim.computational-time-neval",Feel::po::value<int>()->default_value( 0 )," number of evaluation to perform to have the computational time of eim online step" )
        ( "eim.compute-expansion-of-expression",Feel::po::value<bool>()->default_value( false )," use true expression if true, else use the projection of the expression" )
        ;

    return eimoptions;
}
}


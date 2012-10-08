/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file residualestimator.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-07-15
 */
#include <residualestimator.hpp>

#include <feel/feel.hpp>
using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description residualestimatoroptions( "ResidualEstimator options" );
    residualestimatoroptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "dim", po::value<int>()->default_value( 2 ), "dimension of the geometry( 0: all three, 1, 2 or 3" )
    ( "order", po::value<int>()->default_value( 0 ), "order of finite element approximation, 0: execute all registered orders" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 50 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "alpha", Feel::po::value<double>()->default_value( 3 ), "Regularity coefficient for function f" )
    ( "beta", Feel::po::value<double>()->default_value( 1 ), "Coefficient for exponential" )
    ( "fn", Feel::po::value<int>()->default_value( 1 ), "example function to be run" )
    ( "adapt-error-type", Feel::po::value<int>()->default_value( 1 ),"type of error (=1 error indicator, =2 exact error)" )
    ( "adapt-tolerance", Feel::po::value<double>()->default_value( 1e-2 ),"tolerence parameter on the error for mesh adaptation" )
    ( "adapt-hmax", Feel::po::value<double>()->default_value( 2 ),"maximum acceptable h" )
    ( "adapt-hmin", Feel::po::value<double>()->default_value( 1e-5 ),"minimum acceptable h" )
    ( "gmshmodel", Feel::po::value<bool>()->default_value( false ),"enable gmsh model" )
    ( "gmshgeo", Feel::po::value<bool>()->default_value( false ),"enable gmsh model geo file" )

    ;
    return residualestimatoroptions.add( Feel::feel_options() );
}


/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
 
    Application app;


    //app.add( new ResidualEstimator<1,1>() );
    app.add( new ResidualEstimator<2,1>() );
    //app.add( new ResidualEstimator<3,1>() );

    app.run();
}






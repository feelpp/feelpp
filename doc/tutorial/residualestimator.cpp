/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-07-15
 */
#include <residualestimator.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Life::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description residualestimatoroptions("ResidualEstimator options");
    residualestimatoroptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
        ("dim", po::value<int>()->default_value( 0 ), "dimension of the geometry( 0: all three, 1, 2 or 3")
        ("order", po::value<int>()->default_value( 0 ), "order of finite element approximation, 0: execute all registered orders")
        ("shape", Life::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Life::po::value<double>()->default_value( 50 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ("alpha", Life::po::value<double>()->default_value( 3 ), "Regularity coefficient for function f")
        ("fn", Life::po::value<int>()->default_value( 1 ), "example function to be run")
      ("tol", Life::po::value<double>()->default_value(1e-2),"tolerence parameter on the error for mesh adaptation")
        ;
    return residualestimatoroptions.add( Life::life_options() );
}


/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    /**
     * create an application
     */
    /** \code */
    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }
    /** \endcode */

    /**
     * register the simgets
     */
    /** \code */
    app.add( new ResidualEstimator<1,1>( app.vm(), app.about() ) );
    app.add( new ResidualEstimator<2,1>( app.vm(), app.about() ) );
    app.add( new ResidualEstimator<3,1>( app.vm(), app.about() ) );

    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}






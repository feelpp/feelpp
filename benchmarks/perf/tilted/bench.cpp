/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-10-31

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file bench.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-10-31
 */
#include <boost/assign/list_of.hpp>
#include <feel/feelcore/application.hpp>

#include <tilted.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description tiltedoptions( "Tilted options" );
    tiltedoptions.add_options()
#if defined(FEELPP_SOLUTION_1)
        ( "testcase", Feel::po::value<std::string>()->default_value( "solution1" ), "name of the testcase" )
#elif  defined(FEELPP_SOLUTION_2)
        ( "testcase", Feel::po::value<std::string>()->default_value( "solution2" ), "name of the testcase" )
#elif  defined(FEELPP_SOLUTION_3)
        ( "testcase", Feel::po::value<std::string>()->default_value( "solution3" ), "name of the testcase" )
#endif
        ( "faster", Feel::po::value<int>()->default_value( 1 ), "use coupled(0) or default(1) pattern or default/symmetric(2) pattern" )
        ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
        ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0/40 ), "reaction coefficient component" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 400.0 ), "coeff for weak Dirichlet conditions" )
        ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "convection coefficient" )
#if defined(FEELPP_SOLUTION_1) || defined(FEELPP_SOLUTION_2)
        ( "kappa", Feel::po::value<double>()->default_value( 0.1 ), "kappa (testcase 1 or 2 default value)" )
#elif defined(FEELPP_SOLUTION_3)
        ( "kappa", Feel::po::value<double>()->default_value( 30 ), "kappa (testcase 3 default value)" )
#endif
        ( "shear", Feel::po::value<double>()->default_value( 0.0 ), "shear coeff" )
        ( "recombine", Feel::po::value<bool>()->default_value( false ), "recombine triangle into quads" )
        ( "export-matlab", "export matrix and vectors in matlab" )
        ( "no-solve", "dont solve the system" )
        ( "extra-terms", "dont solve the system" )
        ;
    return tiltedoptions.add( Feel::feel_options() )
           .add( Feel::benchmark_options( "2D-CR1-Hypercube" ) ).add( Feel::benchmark_options( "2D-P1-Hypercube" ) ).add( Feel::benchmark_options( "2D-P2-Hypercube" ) );
}


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "tilted" ,
                           "tilted" ,
                           "0.1",
                           "Laplacian equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
namespace Feel
{
extern template class Tilted<2, Lagrange<1, Scalar>, Hypercube>;
extern template class Tilted<2, Lagrange<2, Scalar>, Hypercube>;
extern template class Tilted<2, Lagrange<3, Scalar>, Hypercube>;
extern template class Tilted<2, CrouzeixRaviart<1, Scalar>, Hypercube>;

}

int main( int argc, char** argv )
{

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Environment::changeRepository( boost::format( "%1%/%2%" )
                                   % makeAbout().appName()
                                   % option(_name="testcase").as<std::string>() );

    Application benchmark;

    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
    }

    benchmark.add( new Tilted<2, Lagrange<1, Scalar>, Hypercube>( "2D-P1-Hypercube", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Tilted<2, Lagrange<2, Scalar>, Hypercube>( "2D-P2-Hypercube", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Tilted<2, Lagrange<3, Scalar>, Hypercube>( "2D-P3-Hypercube", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Tilted<2, CrouzeixRaviart<1, Scalar>, Hypercube>( "2D-CR1-Hypercube", benchmark.vm(), benchmark.about() ) );
    benchmark.run();

    //benchmark.printStats( std::cout, boost::assign::list_of( "e.l2" )( "e.h1" )( "n.space" )( "n.matrix" )( "t.init" )( "t.assembly.vector" )( "t.assembly.matrix" )( "t.solver" )( "d.solver" )( "t.integrate" )( "t.export" ) );

    // no h1 yet (need to compute gradient)
    benchmark.printStats( std::cout, boost::assign::list_of( "e.l2" )( "e.h1" )( "n.space" )( "n.matrix" )( "t.init" )( "t.assembly.vector" )( "t.assembly.matrix" )( "t.solver" )( "d.solver" )( "t.integrate" )( "t.export" ), Application::ALL );
}

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-01-15

  Copyright (C) 2013 Université de Strasbourg

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
   \file stokes_pressure.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-01-15
 */

#include <boost/assign/list_of.hpp>
#include <stokes_pressure.hpp>

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
    Feel::po::options_description stokesoptions( "Stokes options" );
    stokesoptions.add_options()
        ( "L", Feel::po::value<double>()->default_value( 1.0 ), "length of the tube" )
        ( "p_in", Feel::po::value<double>()->default_value( 10 ), "pressure at inlet" )
        ( "p_out", Feel::po::value<double>()->default_value( 1 ), "pressure at outlet" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
        ( "geofile", Feel::po::value<std::string>()->default_value( "straighttube.geo" ), "geometry file name" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "bccoefflag", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "eps", Feel::po::value<double>()->default_value( 1e-10 ), "penalisation parameter for lagrange multipliers" )
        ( "gamma-tau", Feel::po::value<double>()->default_value( 10 ), "penalty parameter for tangential velocity" )
        ;
    return stokesoptions.add( Feel::feel_options() ) ;
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
    Feel::AboutData about(
#if !defined( FEELPP_USE_LM )
        "stokes_pressure" ,
        "stokes_pressure" ,
#else
        "stokes_pressure_lm" ,
        "stokes_pressure_lm" ,
#endif
        "0.1",
        "Curvature benchmark using levelset functions",
        Feel::AboutData::License_GPL,
        "Copyright (c) 2012-2013 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
namespace Feel
{
    // 2D
    extern template class Stokes<2, 2, 1>;//P2P1G1
    extern template class Stokes<2, 3, 1>;//P3P2G1
    extern template class Stokes<2, 4, 1>;//P4P3G1
    //3D
    extern template class Stokes<3, 2, 1>;//P2P1G1
    extern template class Stokes<3, 2, 2>;//P2P1G2
    extern template class Stokes<3, 3, 1>;//P3P2G1
    extern template class Stokes<3, 3, 2>;//P3P2G2

}



int main( int argc, char** argv )
{

    using namespace Feel;
    Environment env(_argc=argc,_argv=argv,_desc=makeOptions(),_about=makeAbout());
    Application benchmark;

    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
    }
#if defined(DIM2)
    benchmark.add( new Stokes<2, 2, 1>( "2D-P2P1G1" ) ) ;
    //benchmark.add( new Stokes<2, 3, 1>( "2D-P3P2G1" ) ) ;
#else
    benchmark.add( new Stokes<3, 2, 1>( "3D-P2P1G1" ) ) ;
    //benchmark.add( new Stokes<3, 2, 2>( "3D-P2P1G2" ) ) ;
    //benchmark.add( new Stokes<3, 3, 1>( "3D-P3P2G1" ) ) ;
    //benchmark.add( new Stokes<3, 3, 2>( "3D-P3P2G2" ) ) ;
#endif

    benchmark.setStats( boost::assign::list_of( "e.h1" )( "e.semih1" )( "e.l2" )( "e.output" )("n.space")("d.solver") );

    benchmark.run();
    benchmark.printStats( std::cout );
}



#if 0
int
main( int argc, char** argv )
{

    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(
#if defined( FEELPP_USE_LM )
                     _name="stokes_pressure_lm",
#else
                     _name="stokes_pressure",
#endif
                     _author="Christophe Prud'homme",
                     _email="christophe.prudhomme@feelpp.org") );

    /* define and run application */
#if defined(DIM2)
    Stokes<DIM2,2,1> stokes;
#elif defined(DIM3)
    //Stokes<DIM3,2,1> stokes;
    Stokes<DIM3,3,1> stokes;
#endif
    stokes.run();
}


#endif

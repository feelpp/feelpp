/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-10-31

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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

#include <stokes_kovasnay_curve.hpp>

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
     Feel::po::options_description stokesKovasnayCurveOptions( "Stokes_Kovasnay_Curve options" );
    stokesKovasnayCurveOptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokesKovasnayCurveOptions.add( Feel::feel_options() ) ;
}

/*inline
Feel::po::options_description
makeBenchmarkOptions( std::string const& bench )
{
    Feel::po::options_description benchoptions( "Stokes benchmark options" );
    benchoptions.add_options()
    ( Feel::prefixvm( bench,"bctype" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( Feel::prefixvm( bench,"bccoeff" ).c_str(), Feel::po::value<double>()->default_value( 400.0 ), "coeff for weak Dirichlet conditions" )
    ;
    return benchoptions.add( Feel::benchmark_options( bench ) );
    }*/


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
                          "stokes_kovasnay_curve" ,
                          "stokes_kovasnay_curve" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
namespace Feel
{
    extern template class  Stokes_Kovasnay_Curve<1,1>;
    extern template class  Stokes_Kovasnay_Curve<1,2>;
    extern template class  Stokes_Kovasnay_Curve<1,3>;
    extern template class  Stokes_Kovasnay_Curve<1,4>;

    /* extern template class  Stokes_Kovasnay_Curve<2,1>;
    extern template class  Stokes_Kovasnay_Curve<2,2>;
    extern template class  Stokes_Kovasnay_Curve<2,3>;
    extern template class  Stokes_Kovasnay_Curve<2,4>;

    extern template class  Stokes_Kovasnay_Curve<3,1>;
    extern template class  Stokes_Kovasnay_Curve<3,2>;
    extern template class  Stokes_Kovasnay_Curve<3,3>;
    extern template class  Stokes_Kovasnay_Curve<3,4>;

    extern template class  Stokes_Kovasnay_Curve<4,1>;
    extern template class  Stokes_Kovasnay_Curve<4,2>;
    extern template class  Stokes_Kovasnay_Curve<4,3>;
    extern template class  Stokes_Kovasnay_Curve<4,4>;

    extern template class  Stokes_Kovasnay_Curve<5,1>;
    extern template class  Stokes_Kovasnay_Curve<5,2>;
    extern template class  Stokes_Kovasnay_Curve<5,3>;
    extern template class  Stokes_Kovasnay_Curve<5,4>;*/
}

int main( int argc, char** argv )
{

    /* using namespace Feel;
    Environment env(argc, argv);
    std::ofstream out;
    if ( env.worldComm().rank() == 0 )
        out.open( (boost::format("res-%1%.dat") % env.numberOfProcessors() ).str().c_str() );
    std::vector<std::string> boptions = boost::assign::list_of( "P2P1G1" )( "P2P1G2" )( "P2P1G3" )( "P2P1G4" )
( "P3P2G1" )( "P3P2G2" )( "P3P2G3" )( "P3P2G4" )
( "P4P3G1" )( "P4P3G2" )( "P4P3G3" )( "P4P3G4" )
( "P5P4G1" )( "P5P4G2" )( "P5P4G3" )( "P5P4G4" )
( "P6P5G1" )( "P6P5G2" )( "P6P5G3" )( "P6P5G4" );
    auto cmdoptions = makeOptions();
    BOOST_FOREACH( auto o, boptions )
    {
        cmdoptions.add( makeBenchmarkOptions( o ) );
    }
    Application benchmark;

    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
        }*/
    using namespace Feel;
    using namespace Feel::vf;
    Environment env(_argc=argc,_argv=argv,_desc=makeOptions(),_about=makeAbout());
    std::ofstream out;
    Application benchmark;

    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
    }

    benchmark.add( new Stokes_Kovasnay_Curve<1,1>( "P2P1G1"));
    benchmark.add( new Stokes_Kovasnay_Curve<1,2>( "P2P1G2"));
    /* benchmark.add( new Stokes_Kovasnay_Curve<1,1>( "P2P1G1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<1,2>( "P2P1G2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<1,3>( "P2P1G3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<1,4>( "P2P1G4", benchmark.vm(), benchmark.about() ) );

    benchmark.add( new Stokes_Kovasnay_Curve<2,1>( "P3P2G1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<2,2>( "P3P2G2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<2,3>( "P3P2G3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<2,4>( "P3P2G4", benchmark.vm(), benchmark.about() ) );

    benchmark.add( new Stokes_Kovasnay_Curve<3,1>( "P4P3G1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<3,2>( "P4P3G2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<3,3>( "P4P3G3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<3,4>( "P4P3G4", benchmark.vm(), benchmark.about() ) );

    benchmark.add( new Stokes_Kovasnay_Curve<4,1>( "P5P4G1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<4,2>( "P5P4G2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<4,3>( "P5P4G3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<4,4>( "P5P4G4", benchmark.vm(), benchmark.about() ) );

    benchmark.add( new Stokes_Kovasnay_Curve<5,1>( "P6P5G1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<5,2>( "P6P5G2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<5,3>( "P6P5G3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes_Kovasnay_Curve<5,4>( "P6P5G4", benchmark.vm(), benchmark.about() ) );*/

    benchmark.setStats( boost::assign::list_of( "u.l2" )( "u.h1" ) );
    //std::string "
    benchmark.run();
    benchmark.printStats( std::cout );
    benchmark.printStats( out );
}

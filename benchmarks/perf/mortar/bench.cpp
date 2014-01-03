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

#include <mortar.hpp>

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
    Feel::po::options_description mortaroptions( "Mortar options" );
    mortaroptions.add_options()
    ( "coeff", Feel::po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", Feel::po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),"penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "shear", Feel::po::value<double>()->default_value( 0.0 ), "shear coeff" )
    ( "recombine", Feel::po::value<bool>()->default_value( false ), "recombine triangle into quads" )
    ;
    return mortaroptions.add( Feel::feel_options() )
           .add( Feel::benchmark_options( "2D-P2-P2" ) );//.add( Feel::benchmark_options( "2D-P2-P3" ) ).add( Feel::benchmark_options( "2D-P3-P2" ) )
    // .add( Feel::benchmark_options( "3D-P2-P2" )).add( Feel::benchmark_options( "3D-P3-P3" ));
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "mortar" ,
                           "mortar" ,
                           "0.1",
                           "nD(n=2,3) Mortar using mortar",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@e.ujf-grenoble.fr", "" );
    return about;

}

namespace Feel
{
extern template class MortarBench<2,2,2>;
#if 0
extern template class MortarBench<2,2,3>;
extern template class MortarBench<2,3,2>;
extern template class MortarBench<3,2,2>;
extern template class MortarBench<3,3,3>;
#endif
}

int main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Environment::changeRepository( boost::format( "%1%" )
                                   % makeAbout().appName() );


    Application benchmark;

    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
    }

    // app.add( new MortarBench<2,2,2>( app.vm(), app.about() ) );
    benchmark.add( new MortarBench<2,2,2>( "2D-P2-P2", benchmark.vm(), benchmark.about() ) );
#if  0
    benchmark.add( new MortarBench<2,2,3>( "2D-P2-P3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new MortarBench<2,3,2>( "2D-P3-P2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new MortarBench<3,2,2>( "3D-P2-P2", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new MortarBench<3,3,3>( "3D-P3-P3", benchmark.vm(), benchmark.about() ) );
#endif
    benchmark.run();
    benchmark.printStats( std::cout, boost::assign::list_of( "e.l2" )( "e.h1" )( "t.init" )( "t.assembly" )( "t.solver" ), Application::ALL );
}

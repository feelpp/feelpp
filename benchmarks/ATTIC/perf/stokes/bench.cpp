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

#include <stokes.hpp>

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
        ( "faster", Feel::po::value<int>()->default_value( 2 ), "use coupled(0) or default(1) pattern or default/symmetric(2) pattern" )
        ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
        ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0/40 ), "reaction coefficient component" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 400.0 ), "coeff for weak Dirichlet conditions" )
        ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "convection coefficient" )
        ( "shear", Feel::po::value<double>()->default_value( 0.0 ), "shear coeff" )
        ( "recombine", Feel::po::value<bool>()->default_value( false ), "recombine triangle into quads" )
        ( "testcase", Feel::po::value<std::string>()->default_value( "default" ), "name of the testcase" )
        ( "2D.u_exact_x", Feel::po::value<std::string>()->default_value( "" ), "velocity first component" )
        ( "2D.u_exact_y", Feel::po::value<std::string>()->default_value( "" ), "velocity second component" )
        ( "2D.u_exact_z", Feel::po::value<std::string>()->default_value( "" ), "velocity third component" )
        ( "2D.p_exact", Feel::po::value<std::string>()->default_value( "" ), "" )
        ( "3D.u_exact_x", Feel::po::value<std::string>()->default_value( "" ), "" )
        ( "3D.u_exact_y", Feel::po::value<std::string>()->default_value( "" ), "" )
        ( "3D.u_exact_z", Feel::po::value<std::string>()->default_value( "" ), "" )
        ( "3D.p_exact", Feel::po::value<std::string>()->default_value( "" ), "" )
        ( "export-matlab", "export matrix and vectors in matlab" )
        ( "no-solve", "dont solve the system" )
        ( "extra-terms", "dont solve the system" )
        ;
    return stokesoptions.add( Feel::feel_options() );
}

inline
Feel::po::options_description
makeBenchmarkOptions( std::string const& bench )
{
    Feel::po::options_description benchoptions( "Stokes benchmark options" );
    benchoptions.add_options()
    ( Feel::prefixvm( bench,"bctype" ).c_str(), Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( Feel::prefixvm( bench,"bccoeff" ).c_str(), Feel::po::value<double>()->default_value( 400.0 ), "coeff for weak Dirichlet conditions" )
    ;
    return benchoptions.add( Feel::benchmark_options( bench ) );
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
    Feel::AboutData about("stokes_convergence" ,
                          "stokes_convergence" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2011 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
namespace Feel
{
extern template class Stokes<2, CrouzeixRaviart<1, Vectorial>,Lagrange<0, Scalar,Discontinuous>, Simplex>;
extern template class Stokes<2, CrouzeixRaviart<1, Vectorial>,Lagrange<0, Scalar,Discontinuous>, Hypercube>;
//extern template class Stokes<3, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex>;
extern template class Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Hypercube>;
//extern template class Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Hypercube>;
//extern template class Stokes<3, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Hypercube>;
//extern template class Stokes<3, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Hypercube>;
//extern template class Stokes<3, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>;

extern template class Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Hypercube>;
extern template class Stokes<3, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Hypercube>;
}

int main( int argc, char** argv )
{

    using namespace Feel;
    std::vector<std::string> boptions = boost::assign::list_of( "2D-CR1P0-Simplex" )( "2D-CR1P0-Hypercube" )
        ( "2D-P2P1-Simplex" )( "2D-P2P1-Hypercube" )
        ( "3D-P2P1-Simplex" )( "3D-P2P1-Hypercube" )
        ( "2D-P3P2-Simplex" )( "2D-P3P2-Hypercube" )
        ( "2D-P5P4-Simplex" )( "2D-P5P4-Hypercube" )
        ( "3D-P5P4-Simplex" )( "3D-P5P4-Hypercube" );
    auto cmdoptions = makeOptions();
    BOOST_FOREACH( auto o, boptions )
    {
        cmdoptions.add( makeBenchmarkOptions( o ) );
    }

    Environment env( _argc=argc, _argv=argv,
                     _desc=cmdoptions,
                     _about=makeAbout() );

    Environment::changeRepository( boost::format( "%1%/%2%" )
                                   % makeAbout().appName()
                                   % option(_name="testcase").as<std::string>() );


    std::ofstream out;
    if ( env.worldComm().rank() == 0 )
        out.open( (boost::format("res-%1%.dat") % env.numberOfProcessors() ).str().c_str() );
    Application benchmark;

    //benchmark.add( new Stokes<2, CrouzeixRaviart<1, Vectorial>,Lagrange<0, Scalar,Discontinuous>, Simplex>( "2D-CR1P0-Simplex"  ) );
    //benchmark.add( new Stokes<2, CrouzeixRaviart<1, Vectorial>,Lagrange<0, Scalar,Discontinuous>, Hypercube>( "2D-CR1P0-Hypercube"  ) );
    benchmark.add( new Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>( "2D-P2P1-Simplex" ) );
    //benchmark.add( new Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Hypercube>( "2D-P2P1-Hypercube") );
    benchmark.add( new Stokes<2, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>( "2D-P3P2-Simplex" ) );
    //benchmark.add( new Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Hypercube>( "2D-P2P1-Hypercube") );
    //benchmark.add( new Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>( "2D-P5P4-Simplex" ) );
    //benchmark.add( new Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Hypercube>( "2D-P5P4-Hypercube" ) );

    //benchmark.add( new Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Hypercube>( "3D-P2P1-Hypercube" ) );
    //benchmark.add( new Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>( "3D-P2P1-Simplex") );

    benchmark.setStats( boost::assign::list_of( "e.l2" )( "e.h1" )( "n.space" )( "n.matrix" )( "t.init" )( "t.assembly.rhs" )( "t.assembly.lhs" )( "t.solver" )( "d.solver" ) );
    //benchmark.setStats( boost::assign::list_of( "e.l2" )( "e.h1" )( "n.space" )( "n.matrix" )( "t.init" )( "t.assembly.rhs" )( "t.assembly.lhs" )( "t.solver" ));

    benchmark.run();
    benchmark.printStats( std::cout );
    benchmark.printStats( out );
}

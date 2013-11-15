/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-10-31

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
   \file bench.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-10-31
 */
#include <boost/assign/list_of.hpp>
#include <curvature.hpp>

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
    Feel::po::options_description curvatureoptions( "Curvature options" );
    curvatureoptions.add_options()
    ( "faster", Feel::po::value<int>()->default_value( 1 ), "use coupled(0) or default(1) pattern or default/symmetric(2) pattern" )
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0/40 ), "reaction coefficient component" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 400.0 ), "coeff for weak Dirichlet conditions" )
    ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "convection coefficient" )
    ( "shear", Feel::po::value<double>()->default_value( 0.0 ), "shear coeff" )
    ( "recombine", Feel::po::value<bool>()->default_value( false ), "recombine triangle into quads" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ( "no-solve", "dont solve the system" )
    ( "extra-terms", "dont solve the system" )
    ;
    return curvatureoptions.add( Feel::feel_options() )
        .add(Feel::backend_options("projections"))
        // .add( Feel::benchmark_options( "2D-CR1-Hypercube" ) )
        // .add( Feel::benchmark_options( "2D-P1-Hypercube" ) )
        // .add( Feel::benchmark_options( "2D-P2-Hypercube" ))
        // .add( Feel::benchmark_options( "2D-P1-Simplex" ) )
        // .add( Feel::benchmark_options( "2D-P2-Simplex" ))
        // .add( Feel::benchmark_options( "2D-P3-Simplex" ))
        // .add( Feel::benchmark_options( "2D-P4-Simplex" ))
        // .add( Feel::benchmark_options( "2D-P5-Simplex" ))
    ;
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
    Feel::AboutData about( "curvature" ,
                           "curvature" ,
                           "0.1",
                           "Curvature benchmark using levelset functions",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Vincent Doyeux", "developer", "vincent.doyeux@ujf-grenoble.fr", "" );
    return about;

}
namespace Feel
{
// 2D
extern template class Curvature<2, Lagrange<1, Scalar>, Lagrange<1, Vectorial>, Simplex>;
// extern template class Curvature<2, Lagrange<2, Scalar>, Lagrange<2, Vectorial>, Simplex>;
// extern template class Curvature<2, Lagrange<3, Scalar>, Lagrange<3, Vectorial>, Simplex>;
// extern template class Curvature<2, Lagrange<4, Scalar>, Lagrange<4, Vectorial>, Simplex>;
// extern template class Curvature<2, Lagrange<5, Scalar>, Lagrange<5, Vectorial>, Simplex>;

// extern template class Curvature<2, Lagrange<5, Scalar>, Lagrange<5, Vectorial>, Simplex>;

// 3D
// extern template class Curvature<3, Lagrange<1, Scalar>, Lagrange<1, Vectorial>, Hypercube>;
// extern template class Curvature<3, Lagrange<2, Scalar>, Lagrange<2, Vectorial>, Hypercube>;
// extern template class Curvature<3, Lagrange<1, Scalar>, Lagrange<3, Vectorial>, Simplex>;
// extern template class Curvature<3, Lagrange<2, Scalar>, Lagrange<2, Vectorial>, Simplex>;
// extern template class Curvature<3, Lagrange<3, Scalar>, Lagrange<3, Vectorial>, Simplex>;

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
#if 1
    benchmark.add( new Curvature<2, Lagrange<1, Scalar>, Lagrange<1, Vectorial>, Simplex>( "2D-P1-Simplex") );
    // benchmark.add( new Curvature<2, Lagrange<2, Scalar>, Lagrange<2, Vectorial>, Simplex>( "2D-P2-Simplex") );
    // benchmark.add( new Curvature<2, Lagrange<3, Scalar>, Lagrange<3, Vectorial>, Simplex>( "2D-P3-Simplex") );
    // benchmark.add( new Curvature<2, Lagrange<4, Scalar>, Lagrange<4, Vectorial>, Simplex>( "2D-P4-Simplex" ) );
    // benchmark.add( new Curvature<2, Lagrange<5, Scalar>, Lagrange<5, Vectorial>, Simplex>( "2D-P5-Simplex" ) );
#endif

#if 0
    benchmark.add( new Curvature<2, Lagrange<5, Scalar>, Lagrange<5, Vectorial,Discontinuous>, Simplex>( "2D-P5-Simplex" ) ) ;
#endif
    //    benchmark.add( new Curvature<3, Lagrange<2, Scalar>, Lagrange<2, Vectorial>, Simplex>( "2D-P3-Simplex" ) );

    benchmark.setStats( boost::assign::list_of( "e.nod" )( "e.l2" )( "e.sm" )( "e.hs" )( "e.opt" )("e.ci")("n.space")("n.spacev") );

    benchmark.run();
    benchmark.printStats( std::cout );
}

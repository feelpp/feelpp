/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    Feel::po::options_description stokesoptions("Stokes options");
    stokesoptions.add_options()
        ("penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter")
        ("f", Feel::po::value<double>()->default_value( 0 ), "forcing term")
        ("mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component")
        ("bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions")
        ("export-matlab", "export matrix and vectors in matlab" )
        ("no-solve", "dont solve the system" )
        ;
    return stokesoptions.add( Feel::feel_options() )
        .add( Feel::benchmark_options( "2D-CR1P0" ) ).add( Feel::benchmark_options( "3D-CR1P0" ) )
        .add( Feel::benchmark_options( "2D-P2P1" ) ).add( Feel::benchmark_options( "3D-P2P1" ) )
        .add( Feel::benchmark_options( "2D-P3P2" ) ).add( Feel::benchmark_options( "3D-P3P2" ) )
        .add( Feel::benchmark_options( "2D-P4P3" ) ).add( Feel::benchmark_options( "3D-P4P3" ) )
        .add( Feel::benchmark_options( "2D-P5P4" ) ).add( Feel::benchmark_options( "3D-P5P4" ) );
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
    Feel::AboutData about( "stokes" ,
                           "stokes" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2011 Universite de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
   return about;

}
namespace Feel
{
extern template class Stokes<2, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex>;
extern template class Stokes<3, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex>;
extern template class Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Simplex>;
extern template class Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>;
extern template class Stokes<3, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>;

}

int main( int argc, char** argv )
{

    using namespace Feel;
    Application benchmark( argc, argv, makeAbout(), makeOptions() );
    if ( benchmark.vm().count( "help" ) )
    {
        std::cout << benchmark.optionsDescription() << "\n";
        return 0;
    }
    //benchmark.add( new Stokes<2, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex>( "2D-CR1P0", benchmark.vm(), benchmark.about() ) );
    //benchmark.add( new Stokes<3, CrouzeixRaviart<1, Vectorial,PointSetEquiSpaced>,Lagrange<0, Scalar,Discontinuous>, Simplex>( "3D-CR1P0", benchmark.vm(), benchmark.about() ) );

    //benchmark.add( new Stokes<2, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>( "2D-P2P1", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes<3, Lagrange<2, Vectorial>,Lagrange<1, Scalar>, Simplex>( "3D-P2P1", benchmark.vm(), benchmark.about() ) );
    //benchmark.add( new Stokes<2, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>( "2D-P3P2", benchmark.vm(), benchmark.about() ) );
    //benchmark.add( new Stokes<3, Lagrange<3, Vectorial>,Lagrange<2, Scalar>, Simplex>( "3D-P3P2", benchmark.vm(), benchmark.about() ) );
#if 0
    benchmark.add( new Stokes<2, Lagrange<4, Vectorial>,Lagrange<3, Scalar>, Simplex>( "2D-P4P3", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes<2, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>( "2D-P5P4", benchmark.vm(), benchmark.about() ) );
    benchmark.add( new Stokes<3, Lagrange<5, Vectorial>,Lagrange<4, Scalar>, Simplex>( "3D-P5P4", benchmark.vm(), benchmark.about() ) );
#endif

    benchmark.run();
    benchmark.printStats( std::cout, boost::assign::list_of( "e.l2")("n.space")("t.init")("t.assembly.vector")("t.assembly.matrix" )("t.solver") );
}

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-02-05

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
#include <feel/feelcore/application.hpp>

#include <dar.hpp>

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
    Feel::po::options_description daroptions( "DAR options" );
    daroptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    //        ("g", Feel::po::value<double>()->default_value( 0 ), "boundary term")
    ( "bx", Feel::po::value<double>()->default_value( 1.0 ), "convection X component" )
    ( "by", Feel::po::value<double>()->default_value( 1.0 ), "convection Y component" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "geomap", Feel::po::value<int>()->default_value( 0 ), "type of geomap for integrals" )
    ( "stiff", Feel::po::value<double>()->default_value( 1.0 ), "stiffness parameter of solution" )
    ( "ring", Feel::po::value<bool>()->default_value( 0 ), "0 = square computational domain, 1 = quarter of a ring as computational domain" )
    //("hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrices and vectors in matlab format" )

    ;
    return daroptions.add( Feel::feel_options() );
    //rm.add( Feel::benchmark_options() );
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
    Feel::AboutData about( "dar" ,
                           "dar" ,
                           "0.1",
                           "Diffusion-Advection-Reaction benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;


}
namespace Feel
{
gmsh_ptrtype
createRing( int Dim, int Order, double meshSize, std::string const& convex )
{
    std::ostringstream ostr;
    std::ostringstream nameStr;
    gmsh_ptrtype gmshp( new Gmsh( Dim, Order ) );
    gmshp->setCharacteristicLength( meshSize );
    ostr << gmshp->preamble() << "\n";

    switch ( Dim )
    {
    case 2:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {0.1,0,0,h/2};\n"
             << "Point(2) = {1,0,0,h};\n"
             << "Point(3) = {0,1,0,h};\n"
             << "Point(4) = {0,0.1,0,h/2};\n"
             << "Point(5) = {0,0,0,h/2};\n"
             << "Line(1) = {1,2};\n"
             << "Circle(2) = {2,5,3};\n"
             << "Line(3) = {3,4};\n"
             << "Circle(4) = {4,5,1};\n"
             << "Line Loop(5) = {1,2,3,4};\n"
             << "Plane Surface(6) = {5};\n"
             << "Physical Line(10) = {1};\n"
             << "Physical Line(20) = {2};\n"
             << "Physical Line(30) = {3};\n"
             << "Physical Line(40) = {4};\n"
             //             << "Physical Line(20) = {1,2,4};\n"
             << "Physical Surface(7) = {6};\n";
        nameStr << "ring-" << convex;
        //        fname = __gmsh.generateSquare( "advectiondg2d", meshSize );//
        break;

        // To be added for 3D something like:
        /*    case 3:
                ostr << "h=" << meshSize << ";\n"
                     << "Point(1) = {-1,-1,-1,h};\n"
                     << "Point(2) = {-1, 1,-1,h};\n"
                     << "Point(3) = { 1, 1,-1,h};\n"
                     << "Point(4) = { 1,-1,-1,h};\n"
                     << "Line(1) = {1,4};\n"
                     << "Line(2) = {4,3};\n"
                     << "Line(3) = {3,2};\n"
                     << "Line(4) = {2,1};\n"
                     << "Line Loop(5) = {3,4,1,2};\n"
                     << "Plane Surface(6) = {5};\n"
                     << "Extrude Surface {6, {0,0,2}};\n"
                     << "Physical Surface(10) = {15,23,6,28};\n"
                     << "Physical Surface(20) = {19,27};\n"
                     << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
                     << "Volume(1) = {31};\n"
                     << "Physical Volume(2) = {1};\n";
                nameStr << "cube." << meshSize;
                break;*/
    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}
}

namespace Feel
{
extern template class DAR<1, 1, Continuous, Hypercube>;
extern template class DAR<2, 1, Continuous, Simplex>;
#if 0
extern template class DAR<5, Lagrange<5, Scalar>, Hypercube>;
extern template class DAR<2, Lagrange<1, Scalar>, Simplex>;
extern template class DAR<2, Lagrange<5, Scalar>, Simplex>;
extern template class DAR<3, Lagrange<1, Scalar>, Simplex>;
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

    //benchmark.add( new DAR<1, 1, Continuous, Hypercube>( benchmark.vm(), benchmark.about() ) );
    benchmark.add( new DAR<2, 1, Continuous, Simplex>( benchmark.vm(), benchmark.about() ) );
    benchmark.run();
    benchmark.printStats( std::cout, boost::assign::list_of( "e.l2" )( "n.space" )( "t.init" )( "t.assembly.vector" )( "t.assembly.matrix" )( "t.solver" )( "t.integrate" ), Application::ALL );
}


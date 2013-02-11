/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>
 
  Date 2013-02-11

  Copyright (C) 2008-2013
       Université Joseph Fourier (Grenoble I)
       Université de Strasbourg

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
   \file mymesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>,
                 Guillaume Dollé <guillaume.dolle@math.unistra.fr>
   \date 2013-02-11
 */
#include <feel/feel.hpp>
using namespace Feel;


/**
 * Custom Feel++ options
 */
inline
Feel::po::options_description
makeOptions()
{
    po::options_description myoptions( "MyMesh options" );
    myoptions.add( feel_options() );
    myoptions.add_options()
        ( "hsize",
          po::value<double>()->default_value( 0.1 ),
          "mesh size" )
        ( "shape",
          po::value<std::string>()->default_value( "hypercube" ),
          "shape of the domain (either simplex or hypercube)" )
        ;

    return myoptions;
}


/**
 * Application description
 */
inline
Feel::AboutData
makeAbout()
{
    AboutData about( "mymesh" ,
	             "mymesh" ,
                     "0.2",
                     "Tutorial - Mesh creation",
                      AboutData::License_GPL,
                     "Copyright (c) 2008-2013 Universite Joseph Fourier, Universite de Strasbourg" );

    about.addAuthor( "developers",
                     "Feel++ Consortium",
                     "feelpp-devel@feelpp.org", "" );
    
    return about;
}


/**
 * Entry point
 */
//\code
//# marker_main #
int main( int argc, char** argv )
{
    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );
                     
    // get options
    const int Dim = 2;
    double meshSize = option( _name="hsize" ).as<double>();
    std::string shape = option( _name="shape" ).as<std::string>();

    // choose exec directory
    Environment::changeRepository( boost::format( "doc/manual/tutorial/%1%/%2%-%3%/h_%4%" )
                                   % Environment::about().appName()
                                   % shape
                                   % Dim
                                   % meshSize );


    // create a mesh with GMSH using Feel++ geometry constructor
    auto mesh = createGMSHMesh( _mesh = new Mesh< Simplex<Dim> > ,
                                _desc = domain( _name=( boost::format( "%1%-%2%" ) 
							    % shape % Dim
						       ).str() ,
                                                _shape=shape,
                                                _dim=Dim,
                                                _h=meshSize ) );

    // export results for post processing
    auto e = exporter( _mesh=mesh );
    e->save();

}   // main
//# endmarker_main #
//\endcode

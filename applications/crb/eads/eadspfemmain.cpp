/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file opuseadspfem.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>


#include <eads.hpp>
#include <opusmodelrb.hpp>
#include <feel/feelcrb/pfemapp.hpp>


int
main( int argc, char** argv )
{
    Feel::PFemApp<Feel::OpusModelRB<2,1,2> > app( argc, argv,
            Feel::makeEadsAbout( "eadspfem" ),
            Feel::makeEadsOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    std::vector<double> X( 7 ),Y( 2 );
    X[0]=10; // kic
    X[1]=1e-2; // fluid flow rate
    X[2]=1e6; // Q
    X[3]=100; // c
    X[4]=4e-3; // ea
    //X[4]=3e-3; // ea
    //X[4]=5e-2; // ea
    X[5]=1; // h
    X[6]=2; // polynomial order
    std::cout << "running...\n";
    app.run( X.data(), X.size(), Y.data(), Y.size() );
    //app.run();
}

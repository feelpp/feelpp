/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-02-06

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-02-06
 */
#include <vector>
#include <eads.hpp>
#include <eadsmfemapp.hpp>

int main( int argc, char** argv )
{
    Feel::EadsMFemApp app( argc, argv,Feel::makeEadsAbout( "eadsmfem" ),Feel::makeEadsOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    std::vector<double> X( 7 ),Y( 2 );
    X[0]=app.vm()["ic1.k"].as<double>(); // kic
    X[1]=app.vm()["fluid.flow-rate"].as<double>(); // fluid flow rate
    X[2]=app.vm()["ic1.Q"].as<double>(); // Q
    X[3]=app.vm()["thermal.c"].as<double>(); // c
    X[4]=app.vm()["air.e"].as<double>(); // ea
    X[5]=app.vm()["hsize"].as<double>(); // h
    X[6]=app.vm()["order-temp"].as<int>(); // polynomial order

    app.run( X.data(), X.size(), Y.data(), Y.size() );
}

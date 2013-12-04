/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2009-2013 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file drivencavity.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-05
 */
#include "drivencavity.hpp"

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description  drivencavityoptions( "Driven Cavity problem options" );
     drivencavityoptions.add_options()
         ( "dim", Feel::po::value<int>()->default_value( 2 ), "geometric dimension" )
         ( "continuation", Feel::po::value<bool>()->default_value( 0 ), "use continuation algorithm: =0 (no) =1(yes)" )
         ( "Re", Feel::po::value<double>()->default_value( 100 ), "Reynolds number" )
         ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
         ;
     return  drivencavityoptions.add( Feel::feel_options() ).add( Feel::backend_options( "newtonns" ) );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "driven_cavity" ,
                           "driven_cavity" ,
                           "0.1",
                           "A Steady state incompressible Navier-Stokes solver using a Newton solver",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Ranine Tarabay", "developer", "ranine.a.tarabay@gmail.com", "" );
    return about;

}



int main( int argc, char** argv )
{

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="driven_cavity",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );
    Feel::DrivenCavity<2> DrivenCavity;
    DrivenCavity.run();
}

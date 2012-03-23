/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-04

  Copyright (C) 2009-2012 Universite Joseph Fourier (Grenoble I)

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
   \file convection_main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-04
 */
#include "convection.hpp"


// command line options
inline po::options_description makeOptions()
{
	po::options_description convectionoptions("Convection Options");
	convectionoptions.add_options()
        // Options
        // Format : (nom, type and default value, brief description )
        ("hsize" , po::value<double>()->default_value( 0.1 ) , "mesh size")
        ("fixpointtol", po::value<double>()->default_value( 1e-8 ), "tolerance for the fix point")
        ("gr", po::value<double>()->default_value( 1 ), "nombre de grashof")
        ("pr", po::value<double>()->default_value( 1e-2 ), "nombre de prandtl")
        ("lefttemp", po::value<double>()->default_value( 0.0 ), "temperature on the left side")
        ("newton", "use newton's method")
        ("penalbc",po::value<double>()->default_value( 10.0 ), "penalisation coefficient for the weak boundary conditions")
        ("maxiter_nlin", po::value<int>()->default_value(100), "maximum nonlinearity iteration")
        ("maxiter_solve", po::value<int>()->default_value(100), "maximum solver iteration")
        ("length", po::value<double>()->default_value(1.0), "length of the room")
        ;

	// return the options as well as the feel options
	return convectionoptions.add( feel_options() );
}


// Definition de la fonction qui donne les infos quand l'option --help est passee
inline AboutData
makeAbout()
{
	// Definition de la structure de donnee pour les infos
	AboutData about("Convection",
                    "Convection",
                    "0.1", 				// Version
                    "Natural convection simulation",// Short comment
                    AboutData::License_GPL ,	// Licence
                    "Copyright (c) SQ 2008\nCopyright (c) 2009 Christophe Prud'homme" );// Copyright

	// Informations sur l'auteur
	about.addAuthor("Quinodoz Samuel","Student","samuel.quinodoz@epfl.ch" ,"main developer");
    about.addAuthor("Christophe Prud'homme","Maintainer","christophe.prudhomme@ujf-grenoble.fr" ,"");

	// Retourne les infos
	return about;
}

extern template class Convection<2,1,2>;

int
main( int argc, char** argv )
{
    using namespace Feel;

    const int Order_s(2);
    const int Order_p(1);
    const int Order_t(2);

    /* assertions handling */
    Feel::Assert::setLog( "convection.assert");

    /* define and run application */
    Convection<Order_s, Order_p, Order_t> myconvection( argc, argv, makeAbout(), makeOptions() );
    myconvection.run();
}


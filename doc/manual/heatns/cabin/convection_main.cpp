
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-09-13

  Copyright (C) 2012 Universite de Strasbourg

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
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-09-13
 */
#include "convection.hpp"

// command line options
inline po::options_description makeOptions()
{
    po::options_description convectionoptions( "Convection Options" );
    convectionoptions.add_options()
    // Options
    // Format : (nom, type and default value, brief description )
    ( "output_dir" , po::value<std::string>()->default_value( "cabin" ) , "output directory" )
    ( "adim" , po::value<int>()->default_value( 1 ) , "adimensioned" )
    ( "hsize" , po::value<double>()->default_value( 0.025 ) , "mesh size" )
    ( "fixpointtol", po::value<double>()->default_value( 1e-8 ), "tolerance for the fix point" )
    ( "gr", po::value<double>()->default_value( 1e2 ), "nombre de grashof" )
    ( "rho", po::value<double>()->default_value( 1.177 ),"fluid density" )
    ( "nu", po::value<double>()->default_value( 18.27 ),"kinematic viscosity" )
    ( "k", po::value<double>()->default_value( 2.22e-5 ),"thermal diffusivity" )
    ( "pC", po::value<double>()->default_value( 1100 ),"heat capacity" )
    ( "pr", po::value<double>()->default_value( 1e-2 ), "nombre de prandtl" )
    ( "lefttemp", po::value<double>()->default_value( 0.0 ), "temperature on the left side" )
    ( "newton", "use newton's method" )
    ( "penalbc",po::value<double>()->default_value( 10.0 ), "penalisation coefficient for the weak boundary conditions" )
    ( "weakdir",po::value<int>()->default_value( 1 ),"weak dirichlet" )
    ( "maxiter_nlin", po::value<int>()->default_value( 100 ), "maximum nonlinearity iteration" )
    ( "maxiter_solve", po::value<int>()->default_value( 100 ), "maximum solver iteration" )
    ( "length", po::value<double>()->default_value( 1.0 ), "length of the room" )
    ( "steady",po::value<int>()->default_value( 1 ),"state steady or not" )
    ( "dt",po::value<double>()->default_value( 1e-2 ),"time step" )
    ( "tf",po::value<double>()->default_value( 1 ),"simulation duration" )
    ( "T0",po::value<double>()->default_value( 300 ),"dirichlet condition value" )
    ( "neum",po::value<double>()->default_value( 10 ),"neumann value" );

    // return the options as well as the feel options
    return convectionoptions.add( feel_options() );
};


// Definition de la fonction qui donne les infos quand l'option --help est passee
inline AboutData
makeAbout()
{
    // Definition de la structure de donnee pour les infos
    AboutData about( "Cabin",
                     "Cabin",
                     "0.1", 				// Version
                     "Convection simulation in an airplane cabin",// Short comment
                     AboutData::License_GPL ,	// Licence
                     "Copyright (c) 2012 Christophe Prud'homme" );// Copyright

    // Informations sur l'auteur
    about.addAuthor( "Christophe Prud'homme","Maintainer","prudhomme@unistra.fr" ,"" );

    // Retourne les infos
    return about;
};

int
main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( argc, argv );
    if ( Environment::worldComm().rank() == 0 )
        std::cout << " number of processors : "  << Environment::numberOfProcessors() << "\n";

    const int Order_s( 2 );
    const int Order_p( 1 );
    const int Order_t( 2 );

    /* assertions handling */
    Feel::Assert::setLog( "cabin.assert" );

    /* define and run application */
    Convection myconvection( argc, argv, makeAbout(), makeOptions() );
    myconvection.run();
}
//./cabine_convection --adim==1 --steady=1 --hsize=0.05 --weakdir=0    -snes_max_it 100 -ksp_converged_reason
//-pc_factor_mat_solver_package umfpack -snes_monitor

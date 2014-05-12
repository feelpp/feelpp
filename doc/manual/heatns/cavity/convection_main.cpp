/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-04
 */

#if CRB_SOLVER == 0
    #include "convection.hpp"
#else
    #include "convection_crb.hpp"
    #include <feel/feelcrb/opusapp.hpp>
    #include <feel/feelcrb/crb.hpp>
    #include <feel/feelcrb/crbmodel.hpp>
    #include <feel/feelcrb/crb_trilinear.hpp>
    #include <feel/feelcrb/crbmodeltrilinear.hpp>
//#include <feel/feelcrb/opusapp_heatns.hpp>
#endif

typedef Eigen::VectorXd theta_vector_type;


// command line options
inline po::options_description makeOptions()
{
    po::options_description convectionoptions( "Convection Options" );
    convectionoptions.add_options()
        // Options
        // Format : (nom, type and default value, brief description )
#if (CONVECTION_DIM == 2)
    #if !CRB_SOLVER
        ( "output_dir" , po::value<std::string>()->default_value( "cavity2D" ) , "output directory" )
    #else
        ( "output_dir" , po::value<std::string>()->default_value( "cavity2Dcrb" ) , "output directory" )
    #endif
#else
    #if !CRB_SOLVER
        ( "output_dir" , po::value<std::string>()->default_value( "cavity3D" ) , "output directory" )
    #else
        ( "output_dir" , po::value<std::string>()->default_value( "cavity3Dcrb" ) , "output directory" )
    #endif
#endif
        ( "input_dir" , po::value<std::string>()->default_value( "FEEL/feelopt/doc/manual/heatns/Mesh/" ) , "input directory" )
        ( "readMesh" , po::value<int>()->default_value( 0 ) , "using mesh in file" )
        ( "mesh_name" , po::value<std::string>()->default_value( "domain.msh" ) , "mesh file name" )
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
        ( "use_continuity" , po::value<bool>()->default_value(true), "use continuity method when using Newton" )
        ( "penalbc",po::value<double>()->default_value( 10.0 ), "penalisation coefficient for the weak boundary conditions" )
        ( "weakdir",po::value<int>()->default_value( 1 ),"weak dirichlet" )
        ( "maxiter_nlin", po::value<int>()->default_value( 100 ), "maximum nonlinearity iteration" )
        ( "maxiter_solve", po::value<int>()->default_value( 100 ), "maximum solver iteration" )
        ( "length", po::value<double>()->default_value( 1.0 ), "length of the room" )
        ( "steady",po::value<int>()->default_value( 1 ),"state steady or not" )
        ( "dt",po::value<double>()->default_value( 1e-2 ),"time step" )
        ( "tf",po::value<double>()->default_value( 1 ),"simulation duration" )
        ( "T0",po::value<double>()->default_value( 300 ),"dirichlet condition value" )
        ( "neum",po::value<double>()->default_value( 10 ),"neumann value" )
        ( "Grmin", po::value<double>()->default_value( 10 ), "Grmin" )
        ( "Prmin", po::value<double>()->default_value( 1 ), "Prmin" )
        ( "Grmax", po::value<double>()->default_value( 10 ), "Grmax" )
        ( "Prmax", po::value<double>()->default_value( 1 ), "Prmax" )
        ( "enable-convection-terms",po::value<bool>()->default_value(true), "enable convections terms" )
        ;
//        ( "no-export", "don't export results" );

    // return the options as well as the feel options
    return convectionoptions.add( feel_options() );
}


// Definition de la fonction qui donne les infos quand l'option --help est passee
inline AboutData
makeAbout()
{
    // Definition de la structure de donnee pour les infos
    AboutData about( "Convection",
                     "Convection",
                     "0.1", 				// Version
                     "Natural convection simulation",// Short comment
                     AboutData::License_GPL ,	// Licence
                     "Copyright (c) SQ 2008\nCopyright (c) 2009-2012 Christophe Prud'homme" );// Copyright

    // Informations sur l'auteur
    about.addAuthor( "Quinodoz Samuel","Student","samuel.quinodoz@epfl.ch" ,"main developer" );
    about.addAuthor( "Christophe Prud'homme","Maintainer","christophe.prudhomme@feelpp.org" ,"" );

    // Retourne les infos
    return about;
}


// Definition de la fonction qui donne les infos quand l'option --help est passee
inline AboutData
about()
{
    // Definition de la structure de donnee pour les infos
    AboutData about( "Convection",
                    "Convection",
                    "0.1", 				// Version
                    "Natural convection simulation",// Short comment
                    AboutData::License_GPL ,	// Licence
                    "Copyright (c) SQ 2008\nCopyright (c) 2009-2012 Christophe Prud'homme" );// Copyright

    // Informations sur l'auteur
    about.addAuthor( "Quinodoz Samuel","Student","samuel.quinodoz@epfl.ch" ,"main developer" );
    about.addAuthor( "Christophe Prud'homme","Maintainer","christophe.prudhomme@feelpp.org" ,"" );

    // Retourne les infos
    return about;
}

int
main( int argc, char** argv )
{

    using namespace Feel;

    /* assertions handling */
    Feel::Assert::setLog( "convection.assert" );

    /* define and run application */
#if CRB_SOLVER == 0
    //Feel::Environment env( argc, argv );
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=makeOptions(),
                           _about=makeAbout() );

    if ( Environment::worldComm().rank() == 0 )

        LOG( INFO ) << "CRB_SOLVER = " << CRB_SOLVER ;
        Convection myconvection( argc, argv, makeAbout(), makeOptions() );
        myconvection.run();
#else
        std::cout << "CRB_SOLVER = " << CRB_SOLVER << std::endl;
        //Feel::OpusApp_heatns<ConvectionCrb> myconvectioncrb( argc, argv, makeAbout(), makeOptions()  );
        //Feel::OpusApp<ConvectionCrb> myconvectioncrb( argc, argv, makeAbout(), makeOptions()  );

        Feel::Environment env( _argc=argc, _argv=argv,
                               _desc=opusapp_options("Convection")
                               .add(crbOptions())
                               .add(makeOptions())
                               .add(podOptions())
                               .add(backend_options("backend-primal"))
                               .add(backend_options("backend-dual"))
                               .add(backend_options("backend-l2"))
                               .add(eimOptions())
                               .add(bdf_options("Convection")),
                               _about=makeAbout() );

    if ( Environment::worldComm().rank() == 0 )
        std::cout << " number of processors : "  << Environment::numberOfProcessors() << "\n";

    //Feel::OpusApp<ConvectionCrb , CRB , CRBModel > myconvectioncrb ;
    Feel::OpusApp<ConvectionCrb , CRBTrilinear , CRBModelTrilinear > myconvectioncrb ;
        myconvectioncrb.run();
#endif

}

//./feel_heatns_natural_convection_cavity_2d --config-file=convection.cfg -snes_max_it 100 -ksp_converged_reason -pc_factor_mat_solver_package umfpack -snes_monitor

//./feel_heatns_natural_convection_cavity_2d_crb --config-file=convection.cfg -snes_max_it 100 -ksp_converged_reason -pc_factor_mat_solver_package umfpack -snes_monitor

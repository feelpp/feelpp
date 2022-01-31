/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "../convection.hpp"
#include <feel/feelmor/opusapp.hpp>
#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/crbmodel.hpp>
#include <feel/feelmor/crb_trilinear.hpp>
#include <feel/feelmor/crbmodeltrilinear.hpp>


/**
 * \brief Build the options list for natural convection crb application
 */
inline po::options_description makeOptions()
{
    po::options_description crbcavityoptions( "Convection crb option" );
    crbcavityoptions.add_options()
        ( "adim" , po::value<int>()->default_value( 1 ) , "adimensioned" )
        ( "fixpointtol", po::value<double>()->default_value( 1e-8 ), "tolerance for the fix point" )
        ( "hsize", po::value<double>()->default_value( 0.025 ), "caracterstic mesh size" )
        ( "gr", po::value<double>()->default_value( 1e2 ), "nombre de grashof" )
        ( "rho", po::value<double>()->default_value( 1.25 ),"fluid density" )
        ( "nu", po::value<double>()->default_value( 1.5e-2 ),"kinematic viscosity" )
        ( "k", po::value<double>()->default_value( 2.5e-2 ),"thermal diffusivity" )
        ( "pC", po::value<double>()->default_value( 1000 ),"heat capacity" )
        ( "pr", po::value<double>()->default_value( 0.71 ), "nombre de prandtl" )
        ( "lefttemp", po::value<double>()->default_value( 0.0 ), "temperature on the left side" )
        ( "use_continuity" , po::value<bool>()->default_value(true), "use continuity method when using Newton" )
        ( "penalbc",po::value<double>()->default_value( 10.0 ), "penalisation coefficient for the weak boundary conditions" )
        ( "maxiter_nlin", po::value<int>()->default_value( 100 ), "maximum nonlinearity iteration" )
        ( "maxiter_solve", po::value<int>()->default_value( 100 ), "maximum solver iteration" )
        ( "length", po::value<double>()->default_value( 1.0 ), "length of the room" )
        ( "steady",po::value<int>()->default_value( 1 ),"state steady or not" )
        ( "T0",po::value<double>()->default_value( 1 ),"dirichlet condition value" )
        ( "neum",po::value<double>()->default_value( 1 ),"neumann value" )
        ( "Grmin", po::value<double>()->default_value( 1 ), "Grmin" )
        ( "Prmin", po::value<double>()->default_value( 1 ), "Prmin" )
        ( "Grmax", po::value<double>()->default_value( 10 ), "Grmax" )
        ( "Prmax", po::value<double>()->default_value( 10 ), "Prmax" )
        ;

    return crbcavityoptions.add( feel_options() );
}

/**
 * \brief Short description of the application
 *
 * Can be display with the option '--help'
 */
inline AboutData
makeAbout( std::string const& app_name )
{
    AboutData about( app_name.c_str() );
    return about;
}


int
main( int argc, char** argv )
{
    using namespace Feel;

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options( ConvectionCrb::name() )
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(eimOptions())
                           .add(bdf_options( ConvectionCrb::name() )),
                           _about=makeAbout( ConvectionCrb::name() ) );

    Feel::OpusApp<ConvectionCrb , CRBTrilinear , CRBModelTrilinear > myconvectioncrb ;
    myconvectioncrb.run();
}

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "../convection.hpp"
#include <feel/feelcrb/opusapp.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/crb_trilinear.hpp>
#include <feel/feelcrb/crbmodeltrilinear.hpp>


/**
 * \brief Build the options list for natural convection crb application
 */
inline po::options_description makeOptions()
{
    po::options_description crbcavityoptions( "Convection crb option" );
    crbcavityoptions.add_options()
        ( "rho", po::value<double>()->default_value( 1.17 ),"fluid density" )
        ( "nu", po::value<double>()->default_value( 1.8e-2 ),"kinematic viscosity" )
        ( "k", po::value<double>()->default_value( 2.6e-2 ),"thermal diffusivity" )

        ( "penalbc",po::value<double>()->default_value( 100 ), "penalisation coefficient for the weak boundary conditions" )

        ( "TinletMin", po::value<double>()->default_value( 0.25 ), "mu0 : inlet min temperature" )
        ( "TinletMax", po::value<double>()->default_value( 0.25 ), "mu0 : inlet max temperature" )

        ( "UinletMin", po::value<double>()->default_value( 0.25 ), "mu1 : inlet min velocity" )
        ( "UinletMax", po::value<double>()->default_value( 0.25 ), "mu1 : inlet max velocity" )

        ( "passengers-flux", po::value<double>()->default_value( 35 ), "passengers flux" )

        ( "psiT", po::value<bool>()->default_value( false ), "" )
        ( "delta0", po::value<double>()->default_value( 1 ), "" )
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

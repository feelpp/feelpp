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
        ( "rho", po::value<double>()->default_value( 1.25 ),"fluid density" )
        ( "nu", po::value<double>()->default_value( .5e-2 ),"kinematic viscosity" )
        ( "k", po::value<double>()->default_value( .7e-2 ),"thermal diffusivity" )
        ( "kCpu2", po::value<double>()->default_value( 5 ),"cpu2 thermal diffusivity" )

        ( "use_continuity" , po::value<bool>()->default_value(true), "use continuity method when using Newton" )
        ( "penalbc",po::value<double>()->default_value( 100 ), "penalisation coefficient for the weak boundary conditions" )

        ( "TinletMin", po::value<double>()->default_value( 5 ), "mu0 : min inlet temperature" )
        ( "TinletMax", po::value<double>()->default_value( 20 ), "mu0 : max inlet temperature" )
        ( "Tcpu1Min", po::value<double>()->default_value( 20 ), "mu1 : cpu1 min temperature (dirichlet condtion )" )
        ( "Tcpu1Max", po::value<double>()->default_value( 60 ), "mu1 : cpu1 max temperature (dirichlet conditon)" )
        ( "Tcpu2Min", po::value<double>()->default_value( 20 ), "mu2 : cpu2 min temperature (robin conditon)" )
        ( "Tcpu2Max", po::value<double>()->default_value( 40 ), "mu2 : cpu2 max temperature (robin condition)" )
        ( "UinletMin", po::value<double>()->default_value( 0.25 ), "mu3 : inlet min velocity" )
        ( "UinletMax", po::value<double>()->default_value( 0.25 ), "mu3 : inlet max velocity" )
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

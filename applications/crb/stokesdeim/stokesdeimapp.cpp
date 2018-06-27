/* this file is generated automatically */
#include "stokesdeim.hpp"
#include <feel/feelcrb/opusapp.hpp>
#include <feel/feelcrb/cvgstudy.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description options("stokes deim options");
    options.add_options()
        ( "cvg-study", Feel::po::value<int>()->default_value( 0 ), "");

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("stokesdeim")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(crbSaddlePointOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(options )
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("StokesDeim")),
                           _about=makeStokesDeimAbout( "stokesdeim" ) );


    if ( ioption("cvg-study") )
    {
        Feel::CvgStudy<Feel::StokesDeim,CRBSaddlePoint,CRBModelSaddlePoint> app;
        app.run(ioption("cvg-study"));
    }
    else
    {
        Feel::OpusApp<Feel::StokesDeim,CRBSaddlePoint,CRBModelSaddlePoint> app;
        app.run();
    }

}

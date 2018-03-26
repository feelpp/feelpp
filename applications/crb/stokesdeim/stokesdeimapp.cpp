/* this file is generated automatically */
#include <stokesdeim.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("stokesdeim")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(crbSaddlePointOptions())
                           .add(deimOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("StokesDeim")),
                           _about=makeStokesDeimAbout( "stokesdeim" ) );

    Feel::OpusApp<Feel::StokesDeim,CRBSaddlePoint,CRBModelSaddlePoint> app;
    app.run();
}

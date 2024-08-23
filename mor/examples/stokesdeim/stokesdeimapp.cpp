/* this file is generated automatically */
#include <stokesdeim.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc = argc, _argv = argv,
                           _desc = opusapp_options("stokesdeim")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeStokesDeimOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("StokesDeim")),
                           _about = makeStokesDeimAbout( "stokesdeim" ) );

    Feel::OpusApp<Feel::StokesDeim > app;
    app.run();
}

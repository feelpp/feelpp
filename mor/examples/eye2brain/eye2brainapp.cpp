/* this file is generated automatically */
#include "eye2brain.hpp"
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc = argc, _argv = argv,
                           _desc = opusapp_options("eye2brain")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeEye2BrainOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("Eye2brain")),
                           _about = makeEye2BrainAbout( "eye2brain" ) );

    Feel::OpusApp<Feel::Eye2Brain> app;
    app.run();
}

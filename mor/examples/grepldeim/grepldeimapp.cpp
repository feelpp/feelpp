/* this file is generated automatically */
#include <grepldeim.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc = argc, _argv = argv,
                           _desc = opusapp_options("grepldeim")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeGreplDEIMOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("GreplDEIM")),
                           _about = makeGreplDEIMAbout( "grepldeim" ) );

    Feel::OpusApp<Feel::GreplDEIM<2,2> > app;
    app.run();
}

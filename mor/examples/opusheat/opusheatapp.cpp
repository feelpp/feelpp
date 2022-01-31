/* this file is generated automatically */
#include <opusheat.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("opusheat")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeOpusHeatOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("OpusHeat")),
                           _about=makeOpusHeatAbout( "opusheat" ) );

    Feel::OpusApp<Feel::OpusHeat<false> > app;
    app.run();
}

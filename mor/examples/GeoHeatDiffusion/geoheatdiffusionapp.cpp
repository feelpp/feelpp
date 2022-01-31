/* this file is generated automatically */
#include <geoheatdiffusion.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("geoheatdiffusion")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeGeoHeatDiffusionOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("GeoHeatDiffusion")),
                           _about=makeGeoHeatDiffusionAbout( "geoheatdiffusion" ) );

    Feel::OpusApp<Feel::GeoHeatDiffusion > app;
    app.run();
}

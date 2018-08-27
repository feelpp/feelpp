/* this file is generated automatically */
#include <heatsink2d.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("heatsink2d")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeHeatSink2DOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("HeatSink2D")),
                           _about=makeHeatSink2DAbout( "heatsink2d" ) );

    Feel::OpusApp<Feel::HeatSink2D > app;
    app.run();
}

/* this file is generated automatically */
#include <heat2d.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("heat2d")
                           .add(crbOptions())
                           .add(makeHeat2DOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("Heat2D")),
                           _about=makeHeat2DAbout( "heat2d" ) );

    Feel::OpusApp<Feel::Heat2D > app;
    app.run();
}

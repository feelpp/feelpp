/* this file is generated automatically */
#include <thermalblock-operatorsfree.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("thermalblockfree")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeThermalBlockFreeOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("ThermalBlockFree")),
                           _about=makeThermalBlockFreeAbout( "thermalblockfree" ) );

    Feel::OpusApp<Feel::ThermalBlockFree > app;
    app.run();
}

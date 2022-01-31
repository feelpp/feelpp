/* this file is generated automatically */
#include <thermalblock-minimalversion.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("thermalblockminimalversion")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeThermalBlockMinimalVersionOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("ThermalBlockMinimalVersion")),
                           _about=makeThermalBlockMinimalVersionAbout( "thermalblockminimalversion" ) );

    Feel::OpusApp<Feel::ThermalBlockMinimalVersion > app;
    app.run();
}

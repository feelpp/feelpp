/* this file is generated automatically */
#include <heatshield-minimalversion.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("heatshieldminimalversion3")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeHeatShieldMinimalVersionOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("HeatShieldMinimalVersion")),
                           _about=makeHeatShieldMinimalVersionAbout( "heatshieldminimalversion3" ) );

    Feel::OpusApp<Feel::HeatShieldMinimalVersion<3> > app;
    app.run();
}

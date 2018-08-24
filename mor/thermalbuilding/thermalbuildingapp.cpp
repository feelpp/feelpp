/* this file is generated automatically */
#include <thermalbuilding.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("thermalbuilding")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeThermalBuildingOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("ThermalBuilding")),
                           _about=makeThermalBuildingAbout( "thermalbuilding" ) );

    Feel::OpusApp<Feel::ThermalBuilding > app;
    app.run();
}

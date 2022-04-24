/* this file is generated automatically */
#include <thermalfin.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("thermalfin")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(makeThermalFinOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("ThermalFin")),
                           _about=makeThermalFinAbout( "thermalfin" ) );

    Feel::OpusApp<Feel::ThermalFin > app;
    app.run();
}

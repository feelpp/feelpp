/* this file is generated automatically */
#include <Rbheat.hpp>
#include <feel/feelmor/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::RbHeat> app( argc, argv,
                                     Feel::makeRbHeatAbout( "Rbheat" ),
                                     Feel::makeRbHeatOptions()  );
    app.run();
}

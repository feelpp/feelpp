/* this file is generated automatically */
#include <Cabineheat.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::CabineHeat > app( argc, argv,
                                                      Feel::makeCabineHeatAbout( "Cabineheat" ),
                                                      Feel::makeCabineHeatOptions()  );
    app.run();
}

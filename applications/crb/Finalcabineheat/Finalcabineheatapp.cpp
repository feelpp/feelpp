/* this file is generated automatically */
#include <Finalcabineheat.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::FinalCabineHeat > app( argc, argv,
                                                      Feel::makeFinalCabineHeatAbout( "Finalcabineheat" ),
                                                      Feel::makeFinalCabineHeatOptions()  );
    app.run();
}

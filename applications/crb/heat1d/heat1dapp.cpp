/* this file is generated automatically */
#include <heat1d.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::Heat1D> app( argc, argv,
                                                      Feel::makeHeat1DAbout( "heat1d" ),
                                                      Feel::makeHeat1DOptions()  );
    app.run();
}

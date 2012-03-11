/* this file is generated automatically */
#include <heat1d.hpp>
#include <feel/feelcrb/crbapp.hpp>

int main( int argc, char** argv )
{
    Feel::CrbApp<Feel::Heat1D> app( argc, argv,
                                                      Feel::makeHeat1DAbout( "heat1d" ),
                                                      Feel::makeHeat1DOptions()  );
    app.run();
}

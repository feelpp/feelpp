/* this file is generated automatically */
#include <unsteadyHeat1d.hpp>
#include <feel/feelcrb/crbapp.hpp>

int main( int argc, char** argv )
{
    Feel::CrbApp<Feel::UnsteadyHeat1D> app( argc, argv,
                                                      Feel::makeUnsteadyHeat1DAbout( "unsteadyHeat1d" ),
                                                      Feel::makeUnsteadyHeat1DOptions()  );
    app.run();
}

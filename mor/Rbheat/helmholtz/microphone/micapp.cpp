/* this file is generated automatically */
#include <mic.hpp>
#include <feel/feelcrb/opusapp.hpp>

int main( int argc, char** argv )
{
    Feel::OpusApp<Feel::Microphone> app( argc, argv,
                                         Feel::makeMicrophoneAbout( "mic" ),
                                         Feel::makeMicrophoneOptions()  );
    app.run();
}

#include <cstdlib>

#include <unistd.h>
#include <sstream>

#include "pslogger.hpp"


PsLogger::PsLogger( std::string fileName, std::string format )
    : M_fileName( fileName )
{
    std::stringstream command;
    command << "echo logging output of ps, format: > " << fileName << std::ends;
    system( command.str().c_str() );
    command.str( "" );
    command << "echo " << format << " >> " << fileName << std::ends;
    system( command.str().c_str() );
    command.str( "" );
    command << "ps --no-header -p " << getpid()
            << " -o \"" << format << "\" >> " << M_fileName << std::ends;
    M_command = command.str();
}

void PsLogger::log( std::string logMessage )
{
    if ( logMessage.length() > 0 )
    {
        std::stringstream command;
        command << "echo " << logMessage << " >> " << M_fileName << std::ends;
        system( command.str().c_str() );
    }

    system( M_command.c_str() );
}
